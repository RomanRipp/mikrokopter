#include <iostream>
#include <ros/ros.h>
#include <controls/pidcontroller.h>
#include <mk_msgs/ARMarker.h>
#include <mk_msgs/sensorData.h>
#include <tf/transform_datatypes.h>
#include <geometry_msgs/Pose.h>


int confidence;
geometry_msgs::Pose vispose;
double roll = 0.0;
double pitch = 0.0;
double yaw = 0.0;

void getVisualPos(const mk_msgs::ARMarker msg){
	confidence = msg.confidence;
	vispose = msg.pose.pose;
}

void getSensorOrient(const mk_msgs::sensorData msg){
	//pitch = msg.Nick;
	//roll = msg.Roll;
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
	ros::Publisher motorcomout = nh.advertise<mk_msgs::motorCommands>("/mikrokopter/commands/motor", 5);

	double desiredX = atof(argv[1]);
	double desiredY = atof(argv[2]);
	double desiredZ = atof(argv[3]);

	vispose.position.x = 0;
	vispose.position.y = 0;
	vispose.position.z = 0;
	vispose.orientation.x = 1;
	vispose.orientation.y = 0;
	vispose.orientation.z = 0;
	vispose.orientation.w = 0;

	PIDController pidcontroller;
	//pidcontroller.setDesiredPosition(atof(argv[1]), atof(argv[2]), atof(argv[3]));
	pidcontroller.setDesiredRPY(-1.57);
	pidcontroller.setPIDParameters(atof(argv[4]),atof(argv[5]),atof(argv[6]));
	pidcontroller.setForces(1.0, atof(argv[7]), atof(argv[8]));

	ros::Rate Loop_rate(1000.0/50.0);
	while(ros::ok()){
		if (markerVisible()){
			tf::Quaternion q(vispose.orientation.x, vispose.orientation.y, vispose.orientation.z, vispose.orientation.w);
			tf::Matrix3x3 m(q);
			m.getRPY(roll, pitch, yaw);
			//std::cout<<roll<<" : "<<pitch<<std::endl;
			pidcontroller.setCurrentPosition(vispose.position.x,
												vispose.position.y,
												-vispose.position.z,
												roll,
												pitch,
												yaw);
			pidcontroller.setTime(ros::Time::now().toNSec());
			pidcontroller.setDesiredPosition(desiredX, desiredY, desiredZ);
			//std::cout<<roll<<" : "<<pitch<<std::endl;
			motorcomout.publish(pidcontroller.pid());
		} else {
			//std::cout<<"Marker is not visible"<<std::endl;
			mk_msgs::motorCommands recover;
			recover.roll = 0.0;
			recover.pitch = 0.0;
			recover.yaw = 0.0;
			recover.throttle = 125;
			motorcomout.publish(recover);
		}
		ros::spinOnce();
		Loop_rate.sleep();
	}
	ros::spin();
	return 0;
}
