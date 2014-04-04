#ifndef MK_BEHAVIORS_STATEMACHINE_H
#define MK_BEHAVIORS_STATEMACHINE_H

#include <ros/ros.h>
#include <vector>
#include <std_msgs/Float32.h>
#include <mk_msgs/sensorData.h>
#include <mk_msgs/stateData.h>
#include <mk_msgs/ARMarker.h>
#include <mk_msgs/SetTakeoff.h>
#include <mk_msgs/SetLanding.h>
#include <mk_msgs/SetNavigation.h>
#include <mk_msgs/motorCommands.h>
#include <mk_behaviors/drivercontroller.h>
#include <mk_behaviors/takeoffcontroller.h>
#include <mk_behaviors/navigationcontroller.h>
#include <mk_behaviors/landingcontroller.h>

class StateMachine {

	public:
		/**
		 * Constructor
		 */
		StateMachine();

		/**
		 * Destructor. Close everything nicely.
		 */
		~StateMachine();

		/**
		 * Starts the state machine
		 * This is a blocking method
		 */
		void start(void);

	private:

		enum Behaviors {Takeoff, Navigation, Landing};

		void setupROS(void);
		bool changeBehavior(StateMachine::Behaviors b, double data);
		bool isActive(StateMachine::Behaviors b);
		mk_msgs::SetSerialPoti assignSW(bool careFree, int gps, bool altHld, bool extrernCtrl);

		bool takeoffCallback(
				mk_msgs::SetTakeoff::Request& req, 
				mk_msgs::SetTakeoff::Response& res);
		bool navigationCallback(
				mk_msgs::SetNavigation::Request& req, 
				mk_msgs::SetNavigation::Response& res);
		bool landingCallback(
				mk_msgs::SetLanding::Request& req, 
				mk_msgs::SetLanding::Response& res);
		void getSensorData(
				const mk_msgs::sensorData data);
		void getStateData(
				const mk_msgs::stateData data);
		void getVisualPos(
				const mk_msgs::ARMarker pose);

		ros::NodeHandle nh;
		ros::ServiceServer takeoff_srv, landing_srv, navigation_srv;
		ros::Subscriber sensordata_sub, statedata_sub, vispose_sub;

		ros::Publisher motorcom_pub;
		std::vector<Behaviors> currBehavior;

		DriverController driverControl;
		TakeoffController takeoffControl;
		NavigationController navigationControl;
		LandingController landingControl;
};


#endif /* MK_BEHAVIORS_STATEMACHINE_H */
