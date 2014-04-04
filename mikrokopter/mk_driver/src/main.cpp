#include "driver.h"
#include <ros/ros.h>

int main(int argc, char** argv){
  	ROS_INFO("Hexacopter driver started...");
//	int baudSp = validSpeed(argv[2]);
//	int St = 0;
//	unsigned char buffer=0;

  	ros::init(argc,argv,"wrapper");
  	ros::NodeHandle nh("/mikrokopter/");

	//control structure initialization
  	ros::Subscriber motor_sub = nh.subscribe<mk_msgs::motorCommands>("commands/motor", 1000, motorCallback);
  	ros::ServiceServer setwp_srv = nh.advertiseService("commands/setwp", setwpCallback);
  	ros::ServiceServer listwp_srv = nh.advertiseService("commands/listwp", listwpCallback);
  	ros::ServiceServer setpti_srv = nh.advertiseService("commands/setpoti", setptiCallback);

	//Topic initialization v
	sensordata_pub = nh.advertise<mk_msgs::sensorData>("sensor/data", 1000);
	statedata_pub = nh.advertise<mk_msgs::stateData>("state/data", 1000);

	MkDriver driver = new MkDriver();

	if (driver.openSerial(argv[1], 230400) < 0)
		return -1;

	//Navi Error request
	//encode('e',2,0, 0);
	//ROS_INFO("Request");
	//buffer = NAVI_SEND_INTERVAL; //1 byte sending interval ( in 10ms steps )
	//encode('o',0,&buffer, 1);

	int count = 70;
  	Rate Loop_rate(1000/DRIVER_LOOP_MS);
  	while(ros::ok()){
		count++;
		St = receiveData(); //receiving data from FC
		if (St>0){
			//ROS_INFO("Decode");
			decode(St);
			publish();
		}
		//OSD Data request
		if (count >= 70){
			//ROS_INFO("Request");
			buffer = NAVI_SEND_INTERVAL; //1 byte sending interval ( in 10ms steps )
			encode('o',0,&buffer, 1);
			count = 0;
		}
		spinOnce();
  	}
 	spin();
  	return 0;
}
