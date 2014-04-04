#ifndef MK_BEHAVIORS_PIDController_H
#define MK_BEHAVIORS_PIDController_H

#include <armadillo>
#include "mk_msgs/motorCommands.h"
#include "geometry_msgs/Pose.h"

class PIDController{
	public:
		/*
		/Constructor
		*/
		PIDController();
		/**
		 * Destructor. Close everything nicely.
		 */
		~PIDController();

		/**
		 * The controller
		 */
		mk_msgs::motorCommands pid();
		//controller input:
		arma::vec3 mDesiredPosition;
		arma::vec3 mDesiredVelocity;
		arma::vec3 mDesiredAcceleration;

		//double mGabs;
		double time;
		arma::vec3 mPosition;
		arma::vec3 mPreviousPosition;
		arma::vec3 mRPY;
		arma::vec3 mVelocity;
		double mMass;

		//controller constants:
		double mKD;
		double mKP;
		double mKI;

	private:

		arma::vec3 mErrorInt;
		//controller output:#include <ecl/time.hpp>
		double mF;//Thrust of Quadcopter
		double mRoll;//commanded roll angle
		double mPitch;//commanded roll angle

		//ecl::Stopwatch mClock;
		void limitMinMax(double* value, const double &min, const double &max);
};
#endif /*MK_BEHAVIORS_PIDController_H*/
