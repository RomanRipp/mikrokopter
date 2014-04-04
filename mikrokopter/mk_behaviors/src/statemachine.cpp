#include <mk_behaviors/statemachine.h>
#include <mk_behaviors/params.h>
#include <mk_behaviors/topics.h>
#include <algorithm>

using namespace std;

StateMachine::StateMachine(void) {
	setupROS();
}

StateMachine::~StateMachine(void) {
	takeoffControl.stop();
	nh.shutdown();
}

void StateMachine::start(void) {

	ROS_INFO("Starting behavior state machine");
	while(ros::ok()) {
		
		ros::Duration(LOOP_PERIOD_MS/1000.0).sleep();
		ros::spinOnce();
	}
}

void StateMachine::setupROS(void) {
	takeoff_srv = nh.advertiseService(TAKEOFF_SRV_T, &StateMachine::takeoffCallback, this);
	navigation_srv = nh.advertiseService(NAVIGATION_SRV_T, &StateMachine::navigationCallback, this); 
	landing_srv = nh.advertiseService(LANDING_SRV_T, &StateMachine::landingCallback, this); 
	//Subscribe
	sensordata_sub = nh.subscribe<mk_msgs::sensorData>("/mikrokopter/sensor/data", 1000, &StateMachine::getSensorData, this);
	statedata_sub = nh.subscribe<mk_msgs::stateData>("/mikrokopter/state/data", 1000, &StateMachine::getStateData, this);
	vispose_sub = nh.subscribe<mk_msgs::ARMarker>("/mikrokopter/vision/ar_pose_marker", 1000, &StateMachine::getVisualPos, this);
	//Advertise
	motorcom_pub = nh.advertise<mk_msgs::motorCommands>("/mikrokopter/commands/motor", 5, this);
}

bool StateMachine::takeoffCallback(
		mk_msgs::SetTakeoff::Request& req, 
		mk_msgs::SetTakeoff::Response& res) {

	ROS_INFO("Takeoff request to altitude: %f", req.height);
	if(!changeBehavior(Takeoff, req.height)) {
		ROS_WARN("Could not change behavior to Takeoff");
		res.ack = -1;
	} else {
		res.ack = 1;
	}
	return true;
}

bool StateMachine::navigationCallback(
		mk_msgs::SetNavigation::Request& req, 
		mk_msgs::SetNavigation::Response& res) {

	ROS_INFO("Starting navigation through points 1-%d", req.numofwps);
	if(!changeBehavior(Navigation, req.numofwps)) {
		ROS_WARN("Could not change behavior to Navigation");
		res.ack = -1;
	} else {
		res.ack = 1;
	}
	return true;
}


bool StateMachine::landingCallback(
		mk_msgs::SetLanding::Request& req, 
		mk_msgs::SetLanding::Response& res) {

	ROS_INFO("Landing request to drop down to altitude: %f", req.height);
	landingControl.deflslat = req.lat;
	landingControl.deflslon = req.lon;

	if(!changeBehavior(Landing, req.height)) {
		ROS_WARN("Could not change behavior to Landing");
		res.ack = -1;
	} else {
		res.ack = 1;
	}
	return true;
}

void StateMachine::getSensorData(const mk_msgs::sensorData data){
	takeoffControl.currAlt = data.Altitude;
	landingControl.sensorData = data;
	driverControl.sensorData = data;
	navigationControl.sensorData = data;
}
void StateMachine::getStateData(const mk_msgs::stateData data){
	takeoffControl.position.x = data.x;
	takeoffControl.position.y = data.y;
	takeoffControl.position.z = data.z;
	takeoffControl.orientation.x = data.pitch;
	takeoffControl.orientation.y = data.roll; 
	takeoffControl.orientation.z = data.yaw;
}
void StateMachine::getVisualPos(const mk_msgs::ARMarker pose){
	landingControl.dataind=pose.header.seq;
	landingControl.position.x=pose.pose.pose.position.x;
	landingControl.position.y=pose.pose.pose.position.y;
	landingControl.position.z=pose.pose.pose.position.z;
	landingControl.orientation.x=pose.pose.pose.orientation.x;
	landingControl.orientation.y=pose.pose.pose.orientation.y;
	landingControl.orientation.z=pose.pose.pose.orientation.z;
	landingControl.orientation.w=pose.pose.pose.orientation.w;
}

bool StateMachine::changeBehavior(
		StateMachine::Behaviors b, double data) {

	bool changed = false;
	

	// TODO implement state machine here
	switch(b) {
		case Takeoff:
			if(isActive(Takeoff) || isActive(Landing)) {
				ROS_WARN("Cannot takeoff while in takeoff or landing behavior");
			} else {
				currBehavior.clear();
				currBehavior.push_back(b);
				takeoffControl.motorcomout = motorcom_pub;
				takeoffControl.start(data);
				changed = true;
			}
			break;
				
		case Navigation:
			if(isActive(Takeoff) || isActive(Landing) || isActive(Navigation)) {
				ROS_WARN("Cannot land while in takeoff, landing or navigation behavior");
			} else {
				currBehavior.clear();
				currBehavior.push_back(b);
				navigationControl.start(data);
				changed = true;
			}
			break;


		case Landing:
			if(isActive(Takeoff) || isActive(Landing) || isActive(Navigation)) {
				ROS_WARN("Cannot land while in takeoff, landing or navigation behavior");
			} else {
				currBehavior.clear();
				currBehavior.push_back(b);
				landingControl.motorcomout = motorcom_pub;
				landingControl.start(data);
				changed = true;
			}
			break;
		default:
			ROS_WARN("Incorrect behavior specified.");
			changed = true;
			break;
	}

	return changed;
}

bool StateMachine::isActive(StateMachine::Behaviors b) {

	switch(b){
		case Takeoff: 
			return takeoffControl.running;
		break;
		case Navigation: 
			return navigationControl.running;
		break;
		case Landing: 
			return landingControl.running;
		break;
	}
}
