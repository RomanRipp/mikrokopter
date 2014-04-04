#ifndef MK_BEHAVIORS_CONTROLLER_H
#define MK_BEHAVIORS_CONTROLLER_H

#include <math.h>
#include <mk_msgs/motorCommands.h>

#define RECOVERY_FORWARD_SPEED 0.05
#define RECOVERY_YAW_CHANGE_RATE 0.00001
#define RECOVERY_INITIAL_YAW 0.05

class Controller{
	public:
		/*
		/constructor
		*/
		Controller();
		/*
		/Generate control signal
		*/
		mk_msgs::motorCommands holdposition(double px,   //position of a copter
											double py, 
											double pz,
											double ox, //orientation of a copter
											double oy, 
											double oz, 
											double ow, 
											bool quat, //use quaternions? 
											double setheight); //set height
		/*
		/Recovery if marker is lost
		*/
		mk_msgs::motorCommands recovery(void);
		mk_msgs::motorCommands drop(void);		
		mk_msgs::motorCommands lift(double setheight, double currheight);		

		double PID_P,PID_I,PID_D, zPID_P, zPID_I, zPID_D, TILT_LIMIT;
		int MAX_GAS, MIN_GAS;

	private:
		double err, zerr, perr, pzerr, dz, yawSpeed;
		mk_msgs::motorCommands command;
};
#endif /* MK_BEHAVIORS_CONTROLLER_H */
