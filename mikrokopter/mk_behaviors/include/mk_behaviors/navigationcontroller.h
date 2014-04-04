#ifndef MK_BEHAVIORS_NAVIGATIONCONTROLLER_H
#define MK_BEHAVIORS_NAVIGATIONCONTROLLER_H

#include <mk_behaviors/drivercontroller.h>
#include <mk_msgs/sensorData.h>

class NavigationController {
	public:
		/**
		 * Constructor
		 */
		NavigationController();

		/**
		 * Start takeoff controller
		 */
		void start(double NumOfWPs);

		/**
		 * Stop the takeoff controller
		 */

		mk_msgs::sensorData sensorData;
		void stop(void);
		bool running;
	private:
		void controlloop(void);
		void getrosparams(void);
		
		ros::NodeHandle nh;
		int NumOfWPs;
		DriverController drivercontroller;
};

#endif /* MK_BEHAVIORS_NAVIGATIONCONTROLLER_H */
