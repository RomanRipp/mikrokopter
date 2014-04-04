#include <mk_behaviors/drivercontroller.h>

#define SLEEP_TIME 0.5

//Way point status: 
#define INVALID         0x00
#define NEWDATA         0x01
#define PROCESSED       0x02

DriverController::DriverController(){
	swpoticlt = nh.serviceClient<mk_msgs::SetSerialPoti>("/mikrokopter/commands/setpoti",this);
	setwpclt = nh.serviceClient<mk_msgs::SetGPSWaypoint>("/mikrokopter/commands/setwp",this);
}

DriverController::~DriverController(){
	nh.shutdown();
}


bool DriverController::AssignSW(bool careFree, std::string gps, bool altHld, bool externCtrl){
	mk_msgs::SetSerialPoti poti;	
	if (careFree) poti.request.poti1 = 125;
	else poti.request.poti1 = -125;
	if (altHld) poti.request.poti3 = 125;
	else poti.request.poti3 = -125;
	if (externCtrl) poti.request.poti4 = 125;
	else poti.request.poti4 = -125;
	if (gps == "CH") poti.request.poti2 = 125;
	else if (gps == "PH") poti.request.poti2 = 0;
		else poti.request.poti2 = -125;
	if (swpoticlt.call(poti)){
		//ros::Duration(SLEEP_TIME).sleep();
		return true;
	}else{
		ROS_ERROR("Serial poti switch unsucsessful");
		return false;
	}						
}

bool DriverController::EnableFREE(void){
	if (!this->AssignSW(false, "FREE", false, false)){
		ROS_ERROR("Cannot switch to manual control");
		return false;
	} else {
		ROS_INFO("Manual control: ON");
		return true;
	}
}

bool DriverController::EnableExtCtrl(void){
	bool carefree = !!sensorData.CareFree;
	bool althld = !!sensorData.AltHld;
	if (this->AssignSW(false, "FREE", false, true)){
		ROS_INFO("External control: ON");
		return true;
	} else {
		ROS_ERROR("Cannot switch to External Control, enabling manual control...");
		this->EnableFREE();
		return false;
	} 
}

bool DriverController::EnablePH(void){
	bool swerror = false;
	if (this->AssignSW(false, "FREE", true, false)){
		ROS_INFO("Altitude Hold: ON");
		if (this->AssignSW(false, "PH", true, false)){
			ROS_INFO("Position hold: ON");
			return true;
		} else {swerror = true;}		
	} else {swerror = true;} 

	if(swerror){
		ROS_ERROR("Cannot switch to PH, enabling manual control...");
		this->EnableFREE();
		return false;
	}
}
/* Enable CH with reset
bool DriverController::EnableCH(void){
	bool swerror = false;
	if (this->EnablePH()){
		if (this->AssignSW(true,"PH",true, false)){
			ROS_INFO("Care free: ON");
			if (this->AssignSW(true, "CH", true, false)){
				ROS_INFO("Comming Home: ON");
				return true;			
			} else{swerror = true;}		
		}else{swerror = true;}
	}else{swerror = true;}

	if(swerror){
		ROS_ERROR("Cannot switch to CH, enabling manual control...");
		this->EnableFREE();
		return false;
	}		
}
*/
//Enable CH from PH
bool DriverController::EnableCH(void){
	bool swerror = false;
	if (this->AssignSW(true,"PH",true, false)){
		ROS_INFO("Care free: ON");
		if (this->AssignSW(true, "CH", true, false)){
			ROS_INFO("Comming Home: ON");
			return true;			
		} else{swerror = true;}		
	}else{swerror = true;}

	if(swerror){
		ROS_ERROR("Cannot switch to CH, enabling manual control...");
		this->EnableFREE();
		return false;
	}		
}

bool DriverController::WriteWP(int index, double latitude, double longitude, double altitude, int holdtime, double altituderate, double speed){
	mk_msgs::SetGPSWaypoint point;
	point.request.index = index;
	point.request.latitude = latitude;
	point.request.longitude = longitude;
	point.request.altitude = altitude;
	point.request.status = NEWDATA;
	point.request.holdTime = holdtime;
	point.request.altitudeRate = altituderate;
	point.request.speed = speed;	
	if (setwpclt.call(point)){
//		ros::Duration(SLEEP_TIME).sleep();
		ROS_INFO("Waypoint added");
		return true;
	}else{
		ROS_WARN("Adding waypoint unsucsessful");
		return false;
	}
}

bool DriverController::EraseWPlist(){
	mk_msgs::SetGPSWaypoint eraser;
	eraser.request.index = 0;
	eraser.request.latitude = 0;
	eraser.request.longitude = 0;
	eraser.request.altitude = 0;
	eraser.request.status = INVALID;
	eraser.request.holdTime = 0;
	eraser.request.altitudeRate = 0;
	eraser.request.speed = 0;	
	if (setwpclt.call(eraser)){
		ROS_INFO("Waypoints list erased");
		return true;
	}else{
		ROS_ERROR("Erasing waypoint list unsucsessful");
		return false;
	}						
}
