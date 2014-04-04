#ifndef MK_BEHAVIORS_LANDINGCONTROLLER_H
#define MK_BEHAVIORS_LANDINGCONTROLLER_H

#include <mk_behaviors/controller.h>
#include <mk_behaviors/drivercontroller.h>
#include <mk_msgs/sensorData.h>
#include <mk_msgs/GetRendezvousPoint.h>
#include <EarthCoords.h>

class LandingController {
	public:
		/**
		 * Constructor
		 */
		LandingController();

		/**
		 * Start controller
		 */
		void start(double height);

		/**
		 * Stop the controller
		 */
		void stop(void);

		struct {
			double x,y,z;
		} position;
		struct {
			double x,y,z,w;		
		} orientation;

		mk_msgs::sensorData sensorData;
		ros::Publisher motorcomout;
		
		bool running;
		int dataind;
		double deflslat, deflslon; //Dafault landing spot coordinates
	private:
		void controlloop(void);
		void getrosparams(void);
		bool gpsnavigaion(void);
		bool getLandingcoords(void);
		bool writelswp(void);

		ros::NodeHandle nh;
		double height, lslat, lslon;
		double RECOVER_ALT;
		Controller controller;
		DriverController drivercontroller;
		RSN::EarthCoords earthCoords;
};

#endif /* MK_BEHAVIORS_LANDINGCONTROLLER_H */
