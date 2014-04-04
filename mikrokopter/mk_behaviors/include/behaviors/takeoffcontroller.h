#ifndef MK_BEHAVIORS_TAKEOFFCONTROLLER_H
#define MK_BEHAVIORS_TAKEOFFCONTROLLER_H

#include <mk_behaviors/controller.h>
#include <mk_behaviors/drivercontroller.h>

class TakeoffController {
	public:
		/**
		 * Constructor
		 */
		TakeoffController();

		/**
		 * Start takeoff controller
		 */
		void start(double height);

		/**
		 * Stop the takeoff controller
		 */
		void stop(void);

		struct {
			double x,y,z;
		} position;
		struct {
			double x,y,z,w;		
		} orientation;

		double currAlt;
		bool running;
		ros::Publisher motorcomout;
	private:
		void controlloop(void);
		void getrosparams(void);
		
		ros::NodeHandle nh;
		double height;
		Controller controller;
		DriverController drivercontroller;
};

#endif /* MK_BEHAVIORS_TAKEOFFCONTROLLER_H */
