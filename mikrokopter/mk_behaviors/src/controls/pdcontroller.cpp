#include "controls/pidcontroller.h"

PIDController::PIDController() {
	mTime = 0;

	mPosition.fill(0);
	mPreviousPosition.fill(0);
	mVelocity.fill(0);
	mAngularVelocity.fill(0);
	mRPY.fill(0);
	mPreviousRPY.fill(0);

	mDesiredPosition.fill(0);
	mDesiredRPY.fill(0);
	mDesiredVelocity.fill(0);
	mDesiredAngularVelocity.fill(0);
	mDesiredAcceleration.fill(0);
	mErrorInt.fill(0);
	mAngularErrorInt.fill(0);

	mRoll = 0;
	mPitch = 0;
	motorThrust = 1;
}

PIDController::~PIDController() {

}

void PIDController::setCurrentPosition(double x, double y, double z, double R, double P, double Y){
	mPosition(0)=x;
	mPosition(1)=y;
	mPosition(2)=z;
	mRPY(0)=R;
	mRPY(1)=P;
	mRPY(2)=Y;
}

void PIDController::setDesiredPosition(double x, double y, double z){
	mDesiredPosition(0)=x;
	mDesiredPosition(1)=y;
	mDesiredPosition(2)=z;
}

void PIDController::setDesiredRPY(double yaw){
	mDesiredRPY(0) = 0;
	mDesiredRPY(1) = 0;
	mDesiredRPY(2) = yaw;
}

void PIDController::setTime(double t){
	mTime = t;
}

void PIDController::setForces(double m, double minf, double maxf){
	mMass = m;
	maxF = maxf;
	minF = minf;
	mGabs = 9.8;
}

void PIDController::setPIDParameters(double Kp, double Ki, double Kd){
	mKD = Kd;
	mKP = Kp;
	mKI = Ki;
}

double PIDController::getError(void){
	return arma::norm(mDesiredPosition - mPosition, 2);
}


mk_msgs::motorCommands PIDController::pid(){

	double timediff = (ros::Time::now().toNSec() - mTime)/1000;
	mTime = ros::Time::now().toNSec();
	
	mVelocity = (mPosition - mPreviousPosition) / timediff;
	mAngularVelocity = (mRPY - mPreviousRPY) / timediff;
	mPreviousPosition = mPosition;
	mPreviousRPY = mRPY;

	//############################################################################################
	mErrorInt += (mDesiredPosition - mPosition) * timediff;
	mAngularErrorInt += (mDesiredRPY - mRPY) * timediff;

	err=(pow(px,2)+pow(py,2))/2;
	derr=err-perr;
	v=PID_P*err+PID_I*ierr+PID_D*derr;
	phi = atan2(py,-px);
	if (quat) {
		yaw = atan2(2*(ow*oz+ox*oy),1-2*(pow(oy,2)+pow(oz,2))); 
	} else {
		yaw = oz;
	}
	theta=phi+yaw;
	mPitch=v*sin(theta);
	limitMinMax(&mPitch, -TILT_LIMIT, TILT_LIMIT);	
	//Maximum tilt limiter
	mRoll=v*cos(theta);
	limitMinMax(&mRoll, -TILT_LIMIT, TILT_LIMIT);	
	//PID for altitude controll
	pzerr=zerr;
	zerr=setheight-pz;
	izerr=(zerr+pzerr)/2;
	dzerr=zerr-pzerr;
	dz=dz+zPID_P*zerr+zPID_I*izerr+zPID_D*dzerr;
	if(dz <= MIN_GAS) dz=MIN_GAS;
	else if(dz >= MAX_GAS) dz = MAX_GAS;
	//Publish external control messages 

	mk_msgs::motorCommands command;
	command.roll = mRoll;
	command.pitch = mPitch;
	command.yaw = mYaw;
	command.throttle = mF;
	return command;
}

void PIDController::limitMinMax(double* value, const double& min, const double& max) {
	if (*value < min) {
		*value = min;
	}
	if (*value > max) {
		*value = max;
	}
}
