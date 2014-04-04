#ifndef MK_BEHAVIORS_DRIVERCONTROLLER_H
#define MK_BEHAVIORS_DRIVERCONTROLLER_H

#include <ros/ros.h>
#include <mk_msgs/SetSerialPoti.h>
#include <mk_msgs/sensorData.h>
#include <mk_msgs/SetGPSWaypoint.h>

class DriverController{

	public:
	/*
	//Constructor
	*/
	DriverController();
	/**
	 * Destructor. Close everything nicely.
	 */
	~DriverController();
	/*
	//Serial switch
	*/
	bool AssignSW(bool careFree, std::string gps, bool altHld, bool externCtrl);
	/*
	//Enable manual control
	*/
	bool EnableFREE(void);
	/*
	//Enable external control
	*/
	bool EnableExtCtrl(void);
	/*
	//Enable Position Hold
	*/
	bool EnablePH(void);
	/*
	//Enable Comming Home (GPS navigation)
	*/
	bool EnableCH(void);
	/*
	//Write waypoint
	*/	
	bool WriteWP(int index, double latitude, double longitude, double altitude, int holdtime, double altituderate, double speed);
	/*
	//Erase waypoints
	*/
	bool EraseWPlist();

	mk_msgs::sensorData sensorData;

	private:

	ros::NodeHandle nh;
	ros::ServiceClient swpoticlt;
	ros::ServiceClient setwpclt;
};

#endif /* MK_BEHAVIORS_DRIVERCONTROLLER_H */
