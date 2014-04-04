#include <ros/ros.h>
#include <behaviors/landingcontroller.h>
#include <boost/thread.hpp>

#define LANDING_LOOP_MS 50.0
#define DESCEND_RATE 1.0
#define GROUND_SPEED 2.0

LandingController::LandingController() {
	running = false;
}

void LandingController::start(double h) {

	if(running) {
		ROS_WARN("Already in landing loop");
		return;
	}
	height = h;
	dataind = -1;
	running = true;
	getrosparams();
	if (!gpsnavigaion()){
		ROS_INFO("Landing failed");
		return;	
	}
	ROS_INFO("Strating control loop...");	
	boost::thread runthread(boost::bind(&LandingController::controlloop,this));
}

void LandingController::stop(void) {
	ROS_INFO("Stoppping landing loop");
	nh.shutdown();
	running = false;
}

void LandingController::controlloop(void) {
	bool markervisible=false;
	bool visionlanding=false;
	bool recovery = false;
	bool arrivedtols = false; 
	bool landed=false;
	int preind=dataind;
	
	while(running) {

		//check if pattern is visible
		if(preind==dataind && dataind == -1) {
			markervisible=false;
		} else { 
			markervisible=true;
		}
		preind=dataind;
		dataind=-1;
		if (markervisible){
			if (!visionlanding){ 
				ROS_INFO("Marker detected, enabling vision-based landing...");
				if (drivercontroller.AssignSW(false, "PH", true, true)){
					ROS_INFO("Vision-based landing enabled");
					visionlanding=true;
					recovery=false;
				}else{
					ROS_ERROR("Enabling external control failed, in position hold...");
					return;
				}
			} else {
				motorcomout.publish(controller.holdposition(position.x,
															position.y-0.5, 
															position.z, 
															orientation.x, 
															orientation.y, 
															orientation.z, 
															orientation.w, 
															true, 
															height));
				//TODO determine when to turn motors OFF!
				
			}
		} else {
			//What to do if marker is not visible 
			//Check if arived to lansing spot location
			if(sensorData.CurrentWP == 2 && !!sensorData.TargetReached){
				ROS_INFO("Microcopter arrived to landing spot location");
				arrivedtols = true;
			}else{
				arrivedtols = false;
			}
			if (visionlanding){
				ROS_WARN("Mikrokopter lost sight of marker");
				//ROS_INFO("Starting recovery routine..."); recovery();
				recovery = true;
				visionlanding=false;
				motorcomout.publish(controller.recovery());
			} else if (arrivedtols){
				if (!landed){ 
					ROS_INFO("In position, landing...");
//					if (drivercontroller.EnableExtCtrl()){
						landed=true;
//					}else{
//						return;
//					}
				} else {
				ROS_INFO("Now you can turn the motors OFF");
/*				motorcomout.publish(controller.drop(sensorData.Altitude));
*/
				}							
			}else if (recovery) {
				motorcomout.publish(controller.recovery());			
			}
		}
		ros::Duration(LANDING_LOOP_MS/1000.0).sleep();
	}
}

bool LandingController::gpsnavigaion(void){		
	if (sensorData.PH != 1){
		ROS_INFO("Swithcing to position hold...");
		if (drivercontroller.AssignSW(false, "PH", true, false)){
			ROS_INFO("Position hold enabled");
		}else{
			ROS_ERROR("Failed Switching to PH");
			return false;
		}
	}
	if (writelswp()){
		ROS_INFO("Starting navigation...");
		if (drivercontroller.AssignSW(true, "CH", true, false)){
			ROS_INFO("Navigation started...");
			return true;
		}else{
			ROS_ERROR("Failed enabling navigation");
			return false;
		}
	}else{return false;}
}

bool LandingController::writelswp(void){
	bool failure = false;
	ROS_INFO("Erasing current waypoint list...");
	if (drivercontroller.EraseWPlist() && getLandingcoords()){
		ROS_INFO("Writing landing spot waypoint: Latitude=%f, Longitude=%f",lslat, lslon);
		if(drivercontroller.WriteWP(1, 
									lslat, 
									lslon, 
									height + RECOVER_ALT, 
									0, DESCEND_RATE, 
									GROUND_SPEED) 
		&& drivercontroller.WriteWP(2, 
									lslat, 
									lslon, 
									height, 
									30, 
									DESCEND_RATE, 
									GROUND_SPEED)){
			ROS_INFO("Landing spot way point is written...");
			return true;
		}else{failure = true;}			
	}else{failure = true;}

	if (failure){ 
		ROS_ERROR("Failed sending landing spot way point coordinates");
		return false;
	}
}

bool LandingController::getLandingcoords(void){
	ROS_INFO("Requesting landing spot coordinates...");
	ros::ServiceClient rendezvous_clt = nh.serviceClient<mk_msgs::GetRendezvousPoint>("/rendezvous/point", this);
	mk_msgs::GetRendezvousPoint srv;
	srv.request.currLongitude = sensorData.Longitude;
	srv.request.currLatitude = sensorData.Latitude;

	if(!rendezvous_clt.call(srv)){
		ROS_ERROR("Failed to get coordinates of landing spot, useing default");
		lslat = deflslat;
		lslon = deflslon;
		ROS_INFO("LSLAT: %f < %f",lslat,deflslon);
		return true;	
	}else{
		lslat = srv.response.rendezLatitude;
		lslon = srv.response.rendezLongitude;
		return true;	
	}	
}

void LandingController::getrosparams(void){
	if (!nh.getParam("/mikrokopter/parameter/landing/minGas", controller.MIN_GAS)){
		ROS_WARN("Cannot set min gas parameter, using default value 120");
		controller.MIN_GAS = 120;
	}else{
		ROS_INFO("Minimum gas parameter: %d", controller.MIN_GAS);
	}
	if (!nh.getParam("/mikrokopter/parameter/landing/maxGas", controller.MAX_GAS)){
		ROS_WARN("Cannot set max gas parameter, using default value 150");	
		controller.MAX_GAS = 150;
	}else{
		ROS_INFO("Maximum gas parameter: %d", controller.MAX_GAS);
	}
	if (!nh.getParam("/mikrokopter/parameter/landing/recoverAlt", RECOVER_ALT)){
		ROS_WARN("Cannot set recover altitude parameter, using default value 1");	
		RECOVER_ALT = 6;
	}else{
		ROS_INFO("Recovery altitude parameter: %f", RECOVER_ALT);
	}
	if (!nh.getParam("/mikrokopter/parameter/landing/tiltLimit", controller.TILT_LIMIT)){
		ROS_WARN("Cannot set max tilt parameter, using default value 0.26 rad");	
		controller.TILT_LIMIT = 0.26;
	}else{
		ROS_INFO("Tilt limit parameter: %f", controller.TILT_LIMIT);
	}
	if (!nh.getParam("/mikrokopter/parameter/landing/pidP", controller.PID_P)){
		ROS_WARN("Cannot set PID proprtional component, using default value 50");
		controller.PID_P = 50;	
	}else{
		ROS_INFO("PID proportional parameter: %f", controller.PID_P);
	}
	if (!nh.getParam("/mikrokopter/parameter/landing/pidI", controller.PID_I)){
		ROS_WARN("Cannot set PID integral component, using default value 0");
		controller.PID_I = 0;	
	}else{
		ROS_INFO("PID integral parameter: %f", controller.PID_I);
	}
	if (!nh.getParam("/mikrokopter/parameter/landing/pidD", controller.PID_D)){
		ROS_WARN("Cannot set PID differential component, using default value 0");
		controller.PID_D = 0;	
	}else{
		ROS_INFO("PID differential parameter: %f", controller.PID_D);
	}
	if (!nh.getParam("/mikrokopter/parameter/landing/zpidP", controller.zPID_P)){
		ROS_WARN("Cannot set zPID proprtional component, using default value 50");
		controller.zPID_P = 50;	
	}else{
		ROS_INFO("zPID proportional parameter: %f", controller.zPID_P);
	}
	if (!nh.getParam("/mikrokopter/parameter/landing/zpidI", controller.zPID_I)){
		ROS_WARN("Cannot set zPID integral component, using default value 0");
		controller.zPID_I = 0;	
	}else{
		ROS_INFO("zPID integral parameter: %f", controller.zPID_I);
	}
	if (!nh.getParam("/mikrokopter/parameter/landing/zpidD", controller.zPID_D)){
		ROS_WARN("Cannot set zPID differential component, using default value 0");
		controller.zPID_D = 0;	
	}else{
		ROS_INFO("zPID differential parameter: %f", controller.zPID_D);
	}
}
