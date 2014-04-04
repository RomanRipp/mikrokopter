#include "controls/bkstepcontroller.h"

BkStepController::BkStepController() {
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
	motorThrust = 19.3;
}

BkStepController::~BkStepController() {

}

void BkStepController::setCurrentPosition(double x, double y, double z, double R, double P, double Y){
	mPosition(0)=x;
	mPosition(1)=y;
	mPosition(2)=z;
	mRPY(0)=R;
	mRPY(1)=P;
	mRPY(2)=Y;
}

void BkStepController::setDesiredPosition(double x, double y, double z){
	mDesiredPosition(0)=x;
	mDesiredPosition(1)=y;
	mDesiredPosition(2)=z;
}

void BkStepController::setDesiredRPY(double mRPY(2)){
	mDesiredRPY(0) = 0;
	mDesiredRPY(1) = 0;
	mDesiredRPY(2) = mRPY(2);
}

void BkStepController::setTime(double t){
	mTime = t;
}

void BkStepController::setForces(double m, double f){
	mMass = m;
	maxF = f;
}

void BkStepController::setPIDParameters(double Kp, double Ki, double Kd){
	mKD = Kd;
	mKP = Kp;
	mKI = Ki;
}

double BkStepController::getError(void){
	return arma::norm(mDesiredPosition - mPosition, 2);
}


mk_msgs::motorCommands BkStepController::pid(){

	double timediff = (ros::Time::now().toNSec() - mTime)/1000;
	mTime = ros::Time::now().toNSec();

	mVelocity = (mPosition - mPreviousPosition) / timediff;
	mAngularVelocity = (mRPY - mPreviousRPY) / timediff;
	mPreviousPosition = mPosition;
	mPreviousRPY = mRPY;

	//std::cout<<"Time diff: "<<timediff<<" Velocity: "<<mVelocity<<std::endl;

	//############################################################################################
	mErrorInt += (mDesiredPosition - mPosition) * timediff;
	mAngularErrorInt += (mDesiredRPY - mRPY) * timediff;
/*
	double ctrlX = mDesiredAcceleration(0);
	ctrlX += mKD*(mDesiredVelocity(0) - mVelocity(0));
	ctrlX += mKP*(mDesiredPosition(0) - mPosition(0));
	ctrlX += mKI*mErrorInt(0);

	double ctrlY = mDesiredAcceleration(1);
	ctrlY += mKD*(mDesiredVelocity(1) - mVelocity(1));
	ctrlY += mKP*(mDesiredPosition(1) - mPosition(1));
	ctrlY += mKI*mErrorInt(1);
*/
//	std::cout<<"ctrlX: "<<ctrlX<<" "<<"ctrlY: "<<ctrlY<<std::endl;

	//############################################################################################
	//Thrust
	//############################################################################################
	//mF = -mGabs;
	mF = mDesiredAcceleration(2);
	mF += mKD * (mDesiredVelocity(2) - mVelocity(2));
	mF += mKP * (mDesiredPosition(2) - mPosition(2));
	mF += mKI * mErrorInt(2);
	mF *= mMass / (cos(mRPY(0) * cos(mRPY(1))));

	limitMinMax(&mF, 1, maxF);//maxF N is the maximum force of the QC
//	std::cout<< "Error: " << (mDesiredPosition - mPosition) << std::endl;

//	std::cout << "F: " << mF << " error: " << (mDesiredPosition(2) - mPosition(2)) << std::endl;
	//############################################################################################
	//Roll
	//############################################################################################
	if (mF != 0) {//mF is not allowed to be 0!!!
		//mRoll = -asin(mMass/mF)*((-sin(mRPY(2))*ctrlX + cos(mRPY(2))*ctrlY));
//	} else {
//		mRoll = 0;
		mRoll = 1/(mF*cos(pitch)*cos(roll)*cos(mRPY(2)));
		mRoll *= (-5*y-10*dy-9*mF*sin(pitch)*cos(roll)*cos(mRPY(2));
		mRoll -= 4*mF*dpitch*cos(pitch)*cos(roll)*cos(mRPY(2))+mF*pow(dpitch,2)*sin(pitch)*cos(roll)*cos(mRPY(2));
		mRoll += 2*mF*dmRPY(2)*sin(pitch)*cos(roll)*sin(mRPY(2))+mF*dpitch*dmRPY(2)*cos(pitch)*cos(roll)*sin(mRPY(2));
		mRoll -= mF*dpitch*dmRPY(2)*cos(pitch)*cos(roll)*sin(mRPY(2))-mF*pow(dmRPY(2),2)*sin(pitch)*cos(roll)*cos(mRPY(2)));

	}

	limitMinMax(&mRoll, -M_PI_4/2, M_PI_4/2);//limit to 22.5°
	//std::cout << "mRoll: " << mRoll << " min: " << M_PI_4/2 << std::endl;
//	std::cout << "mRoll: " << mRoll << " error: " << (mDesiredPosition(1) - mPosition(1)) << std::endl;

	//############################################################################################
	//Pitch
	//############################################################################################

	if (mF != 0) {//mF is not allowed to be 0!!!
//		Decimal comPitch = asin( (mass / (comThrust*cos(roll))) * ( cos(mRPY(2))*comX + sin(mRPY(2))*comY));
		mPitch = asin(mMass/(mF*cos(mRPY(0))))*((cos(mRPY(2))*ctrlX + sin(mRPY(2))*ctrlY));
//	} else {
//		mPitch = 0;
		mPitch = 1/(mF*cos(roll)*cos(mRPY(2)))*(-5*x-10*dx-9*mF*sin(roll)*cos(mRPY(2))
			-4*mF*droll*cos(roll)*cos(mRPY(2))+mF*pow(droll,2)*sin(roll)*cos(mRPY(2))
			+2*mF*droll*sin(roll)*sin(mRPY(2))+mF*droll*dmRPY(2)*cos(roll)*sin(mRPY(2))
			-mF*droll*dmRPY(2)*cos(roll)*sin(mRPY(2))-mF*pow(dmRPY(2),2)*sin(roll)*cos(mRPY(2)));
	}

	limitMinMax(&mPitch, -M_PI_4/2, M_PI_4/2);//limit to 22.5°
	//std::cout << "mPitch: " << mPitch << " error: " << (mDesiredPosition(1) - mPosition(1)) << std::endl;

	//############################################################################################
	//Yaw
	//############################################################################################
	if (mF != 0){
		mYaw = mKD*(mDesiredAngularVelocity(2) - mAngularVelocity(2));
		mYaw += mKP*(mDesiredRPY(2) - mRPY(2));
		mYaw += mKI*mAngularErrorInt(2);
	}

	limitMinMax(&mYaw, -M_PI_4/2, M_PI_4/2);//limit to 22.5°
	//std::cout << "Yaw: " << mRPY(2) << std::endl;
	//std::cout << "error: " << mDesiredPosition - mPosition << std::endl;
	std::cout << "abs error: " << arma::norm(mDesiredPosition - mPosition, 2) << std::endl;

	//Backstepping controller
	//double mF = (G+Kp1*(zd-z)+Kd1*(dzd-dz))/(cos(pitch)*cos(roll));

	mk_msgs::motorCommands command;
	command.roll = mRoll;
	command.pitch = mPitch;
	command.yaw = mYaw;
	command.throttle = mF * motorThrust;
	return command;
}

void BkStepController::limitMinMax(double* value, const double& min, const double& max) {
	if (*value < min) {
		*value = min;
	}
	if (*value > max) {
		*value = max;
	}
}
