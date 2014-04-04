#include "mk_controls/pidcontroller.h"

PIDController::PIDController() {
	time = 0;

	mPosition.fill(0);
	mPreviousPosition.fill(0);
	mVelocity.fill(0);
	mRPY.fill(0);

	mDesiredPosition.fill(0);
	mDesiredVelocity.fill(0);
	mDesiredAcceleration.fill(0);
	mErrorInt.fill(0);

	mRoll = 0;
	mPitch = 0;
}

PIDController::~PIDController() {
}

mk_msgs::motorCommands PIDController::pid(){

	double timediff = (ros::Time::now().toNSec() - time)/1000; 
	time = ros::Time::now().toNSec();
	
	//get the current Quadcopter info:
//	mPosition[0] = Position.position.x;
//	mPosition[1] = Position.position.y;
//	mPosition[2] = Position.position.z;
//	mRPY[0] = Position.orientation.x;
//	mVelocity = mQuadcopter->getLinearVelocity();
//	mMass = mQuadcopter->getMass();
	mVelocity = (mPosition - mPreviousPosition) / timediff;
	mPreviousPosition = mPosition;

	//std::cout<<"Time diff: "<<timediff<<" Velocity: "<<mVelocity<<std::endl;

	//############################################################################################
	mErrorInt += (mDesiredPosition - mPosition) * timediff;

	double ctrlX = mDesiredAcceleration(0);
	ctrlX += mKD*(mDesiredVelocity(0) - mVelocity(0));
	ctrlX += mKP*(mDesiredPosition(0) - mPosition(0));
	ctrlX += mKI*mErrorInt(0);

	double ctrlY = mDesiredAcceleration(1);
	ctrlY += mKD*(mDesiredVelocity(1) - mVelocity(1));
	ctrlY += mKP*(mDesiredPosition(1) - mPosition(1));
	ctrlY += mKI*mErrorInt(1);
	
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

	limitMinMax(&mF, 0, 20);//20 N is the maximum force of the QC
//	std::cout<< "Error: " << (mDesiredPosition - mPosition) << std::endl;

//	std::cout << "F: " << mF << " error: " << (mDesiredPosition(2) - mPosition(2)) << std::endl;
	//############################################################################################
	//Roll
	//############################################################################################
	if (mF != 0) {//mF is not allowed to be 0!!!
		mRoll = -asin(mMass/mF)*((-sin(mRPY(2))*ctrlX + cos(mRPY(2))*ctrlY));
//	} else {
//		mRoll = 0;
	}

	limitMinMax(&mRoll, -M_PI_4/2, M_PI_4/2);//limit to 22.5°
//	std::cout << "mRoll: " << mRoll << " error: " << (mDesiredPosition(1) - mPosition(1)) << std::endl;

	//############################################################################################
	//Pitch
	//############################################################################################

	if (mF != 0) {//mF is not allowed to be 0!!!
//		Decimal comPitch = asin( (mass / (comThrust*cos(roll))) * ( cos(yaw)*comX + sin(yaw)*comY));
		mPitch = asin(mMass/(mF*cos(mRPY(0))))*((cos(mRPY(2))*ctrlX + sin(mRPY(2))*ctrlY));
//	} else {
//		mPitch = 0;
	}

	limitMinMax(&mRoll, -M_PI_4/2, M_PI_4/2);//limit to 22.5°
//	std::cout << "mPitch: " << mPitch << " error: " << (mDesiredPosition(1) - mPosition(1)) << std::endl;

	//############################################################################################

//	std::cout << "error: " << mDesiredPosition - mPosition << std::endl;
	std::cout << "abs error: " << arma::norm(mDesiredPosition - mPosition, 2) << std::endl;

	mk_msgs::motorCommands command;
	command.roll = mRoll;
	command.pitch = mPitch;
	command.yaw = 0;
	command.throttle = mF;
	return command;
}

void PIDController::limitMinMax(double* value, const double &min, const double &max){
	if (*value>max){
		*value=max;
	}
	if (*value<min){
		*value=min;
	}
}

