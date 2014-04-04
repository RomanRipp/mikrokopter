#ifndef MK_BEHAVIORS_CONTROLLER_H
#define MK_BEHAVIORS_CONTROLLER_H

#include <math.h>
#include <mk_msgs/motorCommands.h>
#include <geometry_msgs/Pose.h>
#include <tf/transform_datatypes.h>

#define RECOVERY_FORWARD_SPEED 0.05
#define RECOVERY_YAW_CHANGE_RATE 0.00001
#define RECOVERY_INITIAL_YAW 0.05
#define G 127 //Gravitational acceleration

class Controller{
	public:
		/*
		/constructor
		*/
		Controller();
		/*
		/Destructor
		*/
		~Controller();
		/*
		/Generate control signal
		*/
		
		mk_msgs::motorCommands backstepping(geometry_msgs::Pose pose, double zd, double yawd);

		/*
		/Recovery if marker is lost
		*/
		mk_msgs::motorCommands recovery(void);
		mk_msgs::motorCommands drop(void);		
		mk_msgs::motorCommands lift(double setheight, double currheight);		

		geometry_msgs::Pose desiredpose;
		double zPID_P, zPID_I, zPID_D, ACCURACY, RP_COEFF;
		int MAX_GAS, MIN_GAS;

		double Kp1, Kd1, Kd2, Kp2;
		double senspitch, sensroll;

	private:
		double zerr, pzerr, zvel;
		// variables for backtracking:
		double px,py,pz;
		double proll,ppitch,pyaw;
		//Control message
		mk_msgs::motorCommands command;
};
#endif /* MK_BEHAVIORS_CONTROLLER_H */
