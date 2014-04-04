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
	motorThrust = -12.75;
}

PIDController::~PIDController() {

}

void PIDController::setCurrentPosition(double x, double y, double z, double R, double P, double Y){
	mPosition(0)=x;
	mPosition(1)=y;
	mPosition(2)=z;
	//mRPY(0)=R;
	//mRPY(1)=P;
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
	
	if (timediff != 0.0){
		mVelocity = (mPosition - mPreviousPosition) / timediff;
		mAngularVelocity = (mRPY - mPreviousRPY) / timediff;
	}
	mPreviousPosition = mPosition;
	mPreviousRPY = mRPY;

	//std::cout<<"Time diff: "<<timediff<<" Velocity: "<<mVelocity<<std::endl;

	//############################################################################################
	mErrorInt += (mDesiredPosition - mPosition) * timediff;
	mAngularErrorInt += (mDesiredRPY - mRPY) * timediff;

	double ctrlX = mDesiredAcceleration(0);
	ctrlX += mKD*(mDesiredVelocity(0) - mVelocity(0));
	ctrlX += mKP*(mDesiredPosition(0) - mPosition(0));
	ctrlX += mKI*mErrorInt(0);

	double ctrlY = mDesiredAcceleration(1);
	ctrlY += mKD*(mDesiredVelocity(1) - mVelocity(1));
	ctrlY += mKP*(mDesiredPosition(1) - mPosition(1));
	ctrlY += mKI*mErrorInt(1);
	//############################################################################################
	//Thrust
	//############################################################################################
	mF = -mGabs;
	mF += mDesiredAcceleration(2);
	mF += mKD * (mDesiredVelocity(2) - mVelocity(2));
	mF += mKP * (mDesiredPosition(2) - mPosition(2));
	mF += mKI * mErrorInt(2);
	mF *= mMass / (cos(mRPY(0) * cos(mRPY(1))));

	limitMinMax(&mF, minF, maxF);//maxF N is the maximum force of the HC
//	std::cout<< "Error: " << (mDesired2Position - mPosition) << std::endl;

	//std::cout << "F: " << mF << " error: " << (mDesiredPosition(2) - mPosition(2)) << " m: " << mMass / (cos(mRPY(0) * cos(mRPY(1)))) << std::endl;
	//############################################################################################
	//Roll
	//############################################################################################
	if (mF != 0.0) {//mF is not allowed to be 0!!!
		mRoll = -asin(mMass/mF)*((-sin(mRPY(2))*ctrlX + cos(mRPY(2))*ctrlY));
//	} else {
//		mRoll = 0.0;
	}

	limitMinMax(&mRoll, -M_PI_4/2, M_PI_4/2);//limit to 22.5°
	//std::cout << "mRoll: " << mRoll << " : " << asin(mMass/mF) <<" : "<<mMass<<" : "<<mF<< std::endl;
	//std::cout << "mRoll: " << mRoll << " error: " << (mDesiredPosition(1) - mPosition(1)) << std::endl;

	//############################################################################################
	//Pitch
	//############################################################################################

	if (mF != 0.0) {//mF is not allowed to be 0!!!
		mPitch = asin(mMass/(mF*cos(mRPY(0))))*((cos(mRPY(2))*ctrlX + sin(mRPY(2))*ctrlY));
//	} else {
//		mPitch = 0.0;
	}

	limitMinMax(&mPitch, -M_PI_4/2, M_PI_4/2);//limit to 22.5°
	//std::cout << "mPitch: " << mPitch << " error: " << (mDesiredPosition(1) - mPosition(1)) << std::endl;

	//############################################################################################
	//Yaw
	//############################################################################################
	if (mF != 0.0){
		mYaw = mKD*(mDesiredAngularVelocity(2) - mAngularVelocity(2));
		mYaw += mKP*(mDesiredRPY(2) - mRPY(2));
		mYaw += mKI*mAngularErrorInt(2);
	}
	limitMinMax(&mYaw, -2*M_PI_4, 2*M_PI_4);//limit to 45°
	//std::cout << "Yaw: " << mRPY(2) << std::endl;
	//std::cout << "error: " << mDesiredPosition - mPosition << std::endl;
	//std::cout << "abs error: " << arma::norm(mDesiredPosition - mPosition, 2) << std::endl;

	if (mF == 0.0){
		std::cout<<mF<<" = "<<"("<<mGabs
			<<"+"<<mDesiredAcceleration(2)
			<<"+"<<mKD<<"*"<<(mDesiredVelocity(2) - mVelocity(2))
			<<"+"<<mKP<<"*"<<(mDesiredPosition(2) - mPosition(2))
			<<"+"<<mKI<<"*"<<mErrorInt(2)<<")"
			<<"*"<< mMass / (cos(mRPY(0) * cos(mRPY(1))))<<std::endl;
			//exit(11);	
	}

	mk_msgs::motorCommands command;
	command.roll = mRoll;
	command.pitch = - mPitch;
	command.yaw = - mYaw;
	command.throttle = mF * motorThrust;

	//mRPY(0)=mRoll;
	//mRPY(1)= mPitch;
	//mRPY(2)= mYaw;

	return command;
}

void PIDController::limitMinMax(double* value, const double& min, const double& max) {
	if (*value < min) {
		*value = min;
	}
	if (*value > max) {
		*value = max;
	}
	if (isnan(*value)){
		*value = 0.0;
	}
}
