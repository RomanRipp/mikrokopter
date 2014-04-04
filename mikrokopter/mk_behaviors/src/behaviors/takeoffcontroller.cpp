#include <ros/ros.h>
#include <mk_behaviors/takeoffcontroller.h>
#include <boost/thread.hpp>

#define TAKEOFF_LOOP_MS 50.0

TakeoffController::TakeoffController() {
	running = false;
}

void TakeoffController::start(double h) {

	if(running) {
		ROS_WARN("Already in takeoff loop");
		return;
	}

	getrosparams();	
	height = h;
	running = true;
	//Switch to External Control (careFree, gps, altHld, externCtrl)
//	if (drivercontroller.EnableExtCtrl()){
	if (drivercontroller.AssignSW(false, "PH", false, true)){
		ROS_INFO("External control enabled, taking off...");
		boost::thread runthread(boost::bind(&TakeoffController::controlloop,this));
	}else{
		ROS_ERROR("Takeoff failed!");
		return;
	}
}

void TakeoffController::stop(void) {
	//TODO make this a class ModeSwitch
	ROS_INFO("Stoppping takeoff loop.....");	
	if (drivercontroller.AssignSW(false, "PH", true, false)){
		ROS_INFO("Took off, position hold enabled");
		nh.shutdown();
		running = false;
	}else{
		ROS_ERROR("Takeoff failed!");
		return;
	}
}

void TakeoffController::controlloop(void) {
	while(running) {
		motorcomout.publish(controller.lift(height, currAlt));
		if (currAlt >= height) {
			ROS_INFO("Altitude reached");
			stop();
		}
		ros::Duration(TAKEOFF_LOOP_MS/1000.0).sleep();
	}
}

void TakeoffController::getrosparams(void){

	if (!nh.getParam("/mikrokopter/parameter/takeoff/minGas", controller.MIN_GAS)){
		ROS_WARN("Cannot set min gas parameter, using default value 100");
		controller.MIN_GAS = 100;
	} else {
		ROS_INFO("Minimum gas parameter: %d", controller.MIN_GAS);
	}

	if (!nh.getParam("/mikrokopter/parameter/takeoff/maxGas", controller.MAX_GAS)){
		ROS_WARN("Cannot set max gas parameter, using default value 200");	
		controller.MAX_GAS = 200;
	} else {
		ROS_INFO("Maximum gas parameter: %d", controller.MAX_GAS);
	}

	if (!nh.getParam("/mikrokopter/parameter/takeoff/tiltLimit", controller.TILT_LIMIT)){
		ROS_WARN("Cannot set max tilt parameter, using default value 0.26 rad");	
		controller.TILT_LIMIT = 0.26;
	} else {
		ROS_INFO("Tilt limit parameter: %f", controller.TILT_LIMIT);
	}

	if (!nh.getParam("/mikrokopter/parameter/takeoff/pidP", controller.PID_P)){
		ROS_WARN("Cannot set PID proprtional component, using default value 0");
		controller.PID_P = 0;	
	} else {
		ROS_INFO("PID proportional component: %f", controller.PID_P);
	}

	if (!nh.getParam("/mikrokopter/parameter/takeoff/pidI", controller.PID_I)){
		ROS_WARN("Cannot set PID integral component, using default value 0");
		controller.PID_I = 0;	
	} else {
		ROS_INFO("PID integral compoent: %f", controller.PID_I);
	}

	if (!nh.getParam("/mikrokopter/parameter/takeoff/pidD", controller.PID_D)){
		ROS_WARN("Cannot set PID diferential component, using default value 0");
		controller.PID_D = 0;	
	} else {
		ROS_INFO("PID diferential component: %f", controller.PID_D);
	}

	if (!nh.getParam("/mikrokopter/parameter/takeoff/zpidP", controller.zPID_P)){
		ROS_WARN("Cannot set zPID proprtional component, using default value 10");
		controller.zPID_P = 10;	
	} else {
		ROS_INFO("zPID proportional component: %f", controller.zPID_P);
	}

	if (!nh.getParam("/mikrokopter/parameter/takeoff/zpidI", controller.zPID_I)){
		ROS_WARN("Cannot set zPID integral component, using default value 0");
		controller.zPID_I = 0;	
	} else {
		ROS_INFO("zPID integral component: %f", controller.zPID_I);
	}

	if (!nh.getParam("/mikrokopter/parameter/takeoff/zpidD", controller.zPID_D)){
		ROS_WARN("Cannot set zPID differential component, using default value 0");
		controller.zPID_D = 0;		
	} else {
		ROS_INFO("zPID diferential component: %f", controller.zPID_D);
	}
}


