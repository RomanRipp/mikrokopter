#include <mk_behaviors/controller.h>

Controller::Controller(){
	//Controler vatiables initialization
	perr=0; //n-1 step error
	pzerr=0;//n-1 step altitude error
	err=0;	//position error
	zerr=0;	//altitude error
	dz=0;	//horizontal velosity
	yawSpeed = RECOVERY_INITIAL_YAW;
}

mk_msgs::motorCommands Controller::holdposition(double px, 
												double py, 
												double pz,
												double ox,
												double oy, 
												double oz, 
												double ow, 
												bool quat, 
												double setheight){
	double dx=0;
	double dy=0;
	double v=0;
	double ierr=0;
	double derr=0;
	double izerr=0;
	double dzerr=0;
	double yaw=0;
	double phi=0;
	double theta=0;
	mk_msgs::motorCommands msg;
	//PID for going on a landing position
	perr=err;
	err=(pow(px,2)+pow(py,2))/2;
	ierr=(err+perr)/2;
	derr=err-perr;
	v=PID_P*err+PID_I*ierr+PID_D*derr;
	phi = atan2(py,-px);
	if (quat) {
		yaw = atan2(2*(ow*oz+ox*oy),1-2*(pow(oy,2)+pow(oz,2))); 
	} else {
		yaw = oz;
	}
	theta=phi+yaw;
	dy=v*sin(theta);	
	//Maximum tilt limiter
	if(dy <= -TILT_LIMIT) dy=-TILT_LIMIT;
	else if(dy >= TILT_LIMIT) dy = TILT_LIMIT;
	dx=v*cos(theta);
	if(dx <= -TILT_LIMIT) dx=-TILT_LIMIT;
	else if(dx >= TILT_LIMIT) dx = TILT_LIMIT;
	//PID for altitude controll
	pzerr=zerr;
	zerr=setheight-pz;
	izerr=(zerr+pzerr)/2;
	dzerr=zerr-pzerr;
	dz=dz+zPID_P*zerr+zPID_I*izerr+zPID_D*dzerr;
	if(dz <= MIN_GAS) dz=MIN_GAS;
	else if(dz >= MAX_GAS) dz = MAX_GAS;
	//Publish external control messages 
	msg.throttle = dz;
	msg.roll = dx;
	msg.pitch = dy;
	return msg;		
}

mk_msgs::motorCommands Controller::drop(void){
	command.roll = 0;
	command.pitch = 0;
	command.yaw = 0;
	command.throttle = MIN_GAS; //+const maybe? 
	return command;
}

mk_msgs::motorCommands Controller::lift(double setheight, double currheight){
	double izerr=0;
	double dzerr=0;
	pzerr=zerr;
	zerr=setheight-currheight;
	izerr=(zerr+pzerr)/2;
	dzerr=zerr-pzerr;
	dz=dz+zPID_P*zerr+zPID_I*izerr+zPID_D*dzerr;
	if(dz <= MIN_GAS) dz=MIN_GAS;
	else if(dz >= MAX_GAS) dz = MAX_GAS;

	command.roll = 0;
	command.pitch = 0;
	command.yaw = 0;
	command.throttle = dz; 
	return command;
}

mk_msgs::motorCommands Controller::recovery(void){
//	ROS_INFO("angsp %f", yawSpeed);
	yawSpeed = yawSpeed - RECOVERY_YAW_CHANGE_RATE;
	command.pitch = RECOVERY_FORWARD_SPEED;
	command.yaw = yawSpeed;
	if (yawSpeed <= 0) {
		ROS_INFO("Cannot find marker");
		command.pitch = 0;
		command.yaw = 0;
	}
	return command;
}
