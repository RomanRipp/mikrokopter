#include <mk_controls/bkstepcontroller.h>

Controller::Controller(){
	//Controler vatiables initialization
//	perr=0; //n-1 step error
	pzerr=0;//n-1 step altitude error
//	err=0;	//position error
	zerr=0;	//altitude error
	zvel=0;	//horizontal velosity
	//yawSpeed = RECOVERY_INITIAL_YAW;

	px = 0;
	py = 0;
	pz = 0;

	proll = 0;
	ppitch = 0;
	pyaw = 0;
}

Controller::~Controller(){
}

mk_msgs::motorCommands Controller::backstepping(geometry_msgs::Pose pose, double zd, double yawd){

	double x = pose.position.x-desiredpose.position.x;
	double y = pose.position.y-desiredpose.position.y;
	double z = pose.position.z-desiredpose.position.z;
	double dx = x-px;
	double dy = y-py;
	double dz = z-pz;

	double pitch;
	double roll;
	double yaw;	
	tf::Pose tf_pose;
  	tf::poseMsgToTF(pose, tf_pose);
  	tf_pose.getBasis().getRPY(roll, pitch, yaw);
	pitch = - senspitch;
	roll = sensroll;
	//double abs = sqrt(pow(pitch,2));
	//pitch = - pitch/abs*(3.14 - abs);
	double dpitch = pitch - ppitch;
	double droll = roll - proll;
	double dyaw = yaw - pyaw;
	
	
	ROS_INFO("Roll: %f, pitch: %f, yaw: %f, ", roll, pitch, yaw);
	//ROS_INFO("x: %f, dx: %f, y: %f, dy: %f, z: %f, dz: %f", x, dx, y, dy, z, dz);
 	double dzd = 0;
	double dyawd = 0;
//	tf::poseMsgToTF(desiredpose, tf_pose)

	//Backstepping controller
	double u1 = (G+Kp1*(zd-z)+Kd1*(dzd-dz))/(cos(pitch)*cos(roll));
	
	double u3 = 1/(u1*cos(roll)*cos(yaw))*(-5*x-10*dx-9*u1*sin(roll)*cos(yaw)
		-4*u1*droll*cos(roll)*cos(yaw)+u1*pow(droll,2)*sin(roll)*cos(yaw)
		+2*u1*droll*sin(roll)*sin(yaw)+u1*droll*dyaw*cos(roll)*sin(yaw)
		-u1*droll*dyaw*cos(roll)*sin(yaw)-u1*pow(dyaw,2)*sin(roll)*cos(yaw));

	double u2 = 1/(u1*cos(pitch)*cos(roll)*cos(yaw))*(-5*y-10*dy-9*u1*sin(pitch)*cos(roll)*cos(yaw)
		-4*u1*dpitch*cos(pitch)*cos(roll)*cos(yaw)+u1*pow(dpitch,2)*sin(pitch)*cos(roll)*cos(yaw)
		+2*u1*dyaw*sin(pitch)*cos(roll)*sin(yaw)+u1*dpitch*dyaw*cos(pitch)*cos(roll)*sin(yaw)
		-u1*dpitch*dyaw*cos(pitch)*cos(roll)*sin(yaw)-u1*pow(dyaw,2)*sin(pitch)*cos(roll)*cos(yaw));

	double u4 = Kp2*(yawd-yaw)+Kd2*(dyawd-dyaw);
	
	//remember past values
	px = x;	
	py = y;
	pz = z;

	ppitch = pitch;
	proll = roll;
	pyaw = yaw;

	mk_msgs::motorCommands msg;	
	msg.throttle = u1;
	msg.pitch = - RP_COEFF*u2;
	msg.roll = RP_COEFF*u3;
	msg.yaw = - u4;
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
	zvel=zvel+zPID_P*zerr+zPID_I*izerr+zPID_D*dzerr;
	if(zvel <= MIN_GAS) zvel=MIN_GAS;
	else if(zvel >= MAX_GAS) zvel = MAX_GAS;

	command.roll = 0;
	command.pitch = 0;
	command.yaw = 0;
	command.throttle = zvel; 
	return command;
}

mk_msgs::motorCommands Controller::recovery(void){
//	ROS_INFO("angsp %f", yawSpeed);
//	yawSpeed = yawSpeed - RECOVERY_YAW_CHANGE_RATE;
//	command.pitch = RECOVERY_FORWARD_SPEED;
//	command.yaw = yawSpeed;
//	if (yawSpeed <= 0) {
//		ROS_INFO("Cannot find marker");
//		command.pitch = 0;
//		command.yaw = 0;
//	}
	command.yaw = 0;
	command.pitch = 0; 
	command.roll = 0;
	command.throttle = 100;
	return command;
}
