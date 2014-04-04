#include <iostream>
#include <ros/ros.h>
#include <controls/pidcontroller.h>
#include <mk_msgs/ARMarker.h>
#include <mk_msgs/sensorData.h>
#include <tf/transform_datatypes.h>
#include <geometry_msgs/Pose.h>

int confidence;
geometry_msgs::Pose vispose;
double roll;
double pitch;
double desiredX;
double desiredY;

void getVisualPos(const mk_msgs::ARMarker msg){
	confidence = msg.confidence;
	vispose = msg.pose.pose;
}

void getSensorOrient(const mk_msgs::sensorData msg){
	pitch = msg.Nick;
	roll = msg.Roll;
}

void getUGVPose(const geometry_msgs::Point msg){
	desiredX = msg.x;
	desiredY = msg.y;
}

bool markerVisible(){
	if (confidence == 0) {
		return false;
	} else {
		return true;
	}
}

int main(int argc, char** argv) {
	ros::init(argc,argv,"mk_behaviors");
	ros::NodeHandle nh;
	ros::Subscriber vispose_sub = nh.subscribe<mk_msgs::ARMarker>("/mikrokopter/vision/ar_pose_marker", 1000, &getVisualPos);
	ros::Subscriber sensororient = nh.subscribe<mk_msgs::sensorData>("/mikrokopter/sensor/data",1000,&getSensorOrient);
	ros::Subscriber ugvposition = nh.subscribe<geometry_msgs::Point>("/mikrokopter/desired",1000,&getUGVPose);
	ros::Publisher motorcomout = nh.advertise<mk_msgs::motorCommands>("/mikrokopter/commands/motor", 5);


	double desiredZ = atof(argv[3]);

	vispose.position.x = 0;
	vispose.position.y = 0;
	vispose.position.z = 0;
	vispose.orientation.x = 1;
	vispose.orientation.y = 0;
	vispose.orientation.z = 0;
	vispose.orientation.w = 0;

	PIDController pidcontroller;
	pidcontroller.setDesiredRPY(-1.57);
	pidcontroller.setPIDParameters(atof(argv[4]),atof(argv[5]),atof(argv[6]));
	pidcontroller.setForces(2.2, atof(argv[7]), atof(argv[8]));

	while(ros::ok()){
		pidcontroller.setCurrentPosition(vispose.position.x,
											vispose.position.y,
											vispose.position.z,
											roll,
											pitch,
											tf::getYaw(vispose.orientation));
		pidcontroller.setTime(ros::Time::now().toNSec());

		double DesiredPosModRadius = 4;
		double DesiredPosModFreq = 0.1;
		double DesiredPosModStartTime = 0;
		desiredX += 0.1; -DesiredPosModRadius*cos(2*M_PI*DesiredPosModFreq*(DesiredPosModStartTime - ros::Time::now().toSec())) + DesiredPosModRadius;
		desiredY = 1;-DesiredPosModRadius*sin(2*M_PI*DesiredPosModFreq*(DesiredPosModStartTime - ros::Time::now().toSec()));

		desiredZ = - 10;
		pidcontroller.setDesiredPosition(desiredX, desiredY, desiredZ);
		motorcomout.publish(pidcontroller.pid());
		ros::Duration(0.05).sleep();
		ros::spinOnce();
	}

	return 0;
}
