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
		 * Setting the controller parameters
		 */
		void setCurrentPosition(double x, double y, double z, double R, double P, double Y);
		/**
		 * setting the desired position
		 */
		void setDesiredPosition(double x, double y, double z);
		/**
		 * ser desired yaw
		 */
		void setDesiredRPY(double yaw);
		/**
		 * setting the mass of the UAV
		 */
		void setForces(double m, double minf, double maxf);
		/**
		 * setting the parameters of pid controller
		 */
		void setPIDParameters(double Kp, double Ki, double Kd);
		/**
		 * the time is needed for velocity estimation
		 */
		void setTime(double t);
		/**
		 * return error value
		 */
		double getError(void);
		/**
		 * The controller
		 */
		mk_msgs::motorCommands pid();

	private:

		//controller input:
		arma::vec3 mDesiredPosition;
		arma::vec3 mDesiredRPY;
		arma::vec3 mDesiredVelocity;
		arma::vec3 mDesiredAngularVelocity;
		arma::vec3 mDesiredAcceleration;

		//double mGabs;
		double mTime;
		arma::vec3 mPosition;
		arma::vec3 mPreviousPosition;
		arma::vec3 mRPY;
		arma::vec3 mPreviousRPY;
		arma::vec3 mVelocity;
		arma::vec3 mAngularVelocity;
		double mMass;
		double mGabs;
		double maxF;
		double minF;
		double motorThrust;

		//controller constants:
		double mKD;
		double mKP;
		double mKI;

		arma::vec3 mErrorInt;
		arma::vec3 mAngularErrorInt;
		double mF;//Thrust of Quadcopter
		double mRoll;//commanded roll angle
		double mPitch;//commanded roll angle
		double mYaw;//commanded yaw angle

		//ecl::Stopwatch mClock;
		void limitMinMax(double* value, const double &min, const double &max);
};
#endif /*MK_BEHAVIORS_PIDController_H*/
