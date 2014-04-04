#include <ros/ros.h>
#include <mk_behaviors/navigationcontroller.h>
#include <boost/thread.hpp>

#define NAVIGATION_LOOP_MS 50.0

NavigationController::NavigationController() {
	running = false;
}

void NavigationController::start(double n) {

	if(running) {
		ROS_WARN("Already in navigation loop");
		return;
	}
	NumOfWPs = n; 
	getrosparams();	
	if (drivercontroller.AssignSW(true, "CH", true, false)){
		ROS_INFO("Navigation started...");
		running = true;
		boost::thread runthread(boost::bind(&NavigationController::controlloop,this));
	}else{
		ROS_ERROR("Failed enabling navigation");
		return;
	}
}

void NavigationController::stop(void) {
	ROS_INFO("Stoppping naviagtion loop.....");	
	running = false;
	if (drivercontroller.AssignSW(false, "PH", true, false)){
		ROS_INFO("Navigation complete, position hold enabled");
	}else{
		ROS_ERROR("Failed enabling position hold");
	}
	nh.shutdown();
	return;
}

void NavigationController::controlloop(void) {
	while(running) {
	//TODO Navigation control here
	if (sensorData.CurrentWP >= NumOfWPs && !!sensorData.TargetReached){
		stop();
		}
	}
}

void NavigationController::getrosparams(void){
/*
	if (!nh.getParam("/mikrokopter/parameter/", PARAMETER)){
		ROS_WARN("Cannot set parameter, using default value ");
		PARAMETER = 100;
	} else {
		ROS_INFO("Parameter: %d", PARAMETER);
	}
*/
}


