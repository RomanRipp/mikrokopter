#include <ros/ros.h>
#include <mk_msgs/motorCommands.h>
#include <mk_msgs/naviCommands.h>
#include <iostream>
#include <stdio.h>
#include <termios.h>
#include <unistd.h>
#include <sys/select.h>
#include <stropts.h>
#include <asm/ioctls.h>

using namespace std;
using namespace ros;

char c;

int mykbhit() {
    static const int STDIN = 0;
    static bool initialized = false;

    if (! initialized) {
        // Use termios to turn off line buffering
        termios term;
        tcgetattr(STDIN, &term);
        term.c_lflag &= ~ICANON;
        tcsetattr(STDIN, TCSANOW, &term);
        setbuf(stdin, NULL);
        initialized = true;
    }

    int bytesWaiting;
    ioctl(STDIN, FIONREAD, &bytesWaiting);
    return bytesWaiting;
}


int mygetch() {
  struct termios oldt,
                 newt;
  int            ch;
  tcgetattr( STDIN_FILENO, &oldt );
  newt = oldt;
  newt.c_lflag &= ( ICANON | ECHO );
  tcsetattr( STDIN_FILENO, TCSANOW, &newt );
  ch = getchar();
  tcsetattr( STDIN_FILENO, TCSANOW, &oldt );
  return ch;
}


int main(int argc, char** argv){
	int i;	
	ROS_INFO("mikrokopter control node");
	init(argc,argv,"TelOp");
	NodeHandle nh;
	Publisher motorcomout = nh.advertise<mk_msgs::motorCommands>("/mikrokopter/commands/motor",5);
	Publisher navicomout = nh.advertise<mk_msgs::naviCommands>("/mikrokopter/commands/navi",5);
	Rate Loop_rate(10);
	int Gas = 126;
	float Roll = 0;
	float Pitch = 0;
	float Yaw = 0;
	double STEP = 0.01;
	int latitude = atoi(argv[2]);
	int longitude = atoi(argv[1]);
	int altitude = atoi(argv[3]);
	int mode = atoi(argv[4]);
	int index = 1;

	sleep(1);

	//Wright waypoints: 
	mk_msgs::naviCommands msg;
	msg.status = 1;
	msg.heading = 1;
	msg.toleranceRadius = 1;
	msg.holdTime = 1;
	msg.eventFlag = 0;
	msg.type = 0;
	msg.eventChannel = 1;
	msg.altitudeRate = 1;
	msg.speed = 1;
	msg.camAngle = 1;

	//WP#1:
	msg.longitude = longitude;
	msg.latitude = latitude;
	msg.altitude = altitude;		
	msg.index = index;

	cout<<"Latitude: "<<msg.latitude<<" Longitude: "<<msg.longitude<<" Altitude: "<<msg.altitude<<" WP"<<msg.index<<endl;			
	navicomout.publish(msg);
	
	//WP2:
	index++;
	msg.longitude = -93.2329895;
	msg.latitude = 44.9746007;
	msg.altitude = 10;		
	msg.index = index;
	cout<<"Latitude: "<<msg.latitude<<" Longitude: "<<msg.longitude<<" Altitude: "<<msg.altitude<<" WP"<<msg.index<<endl;			
	navicomout.publish(msg);

	//WP3:
	index++;
	msg.longitude = -93.2329895;
	msg.latitude = 44.9746000;
	msg.altitude = 10;		
	msg.index = index;
	cout<<"Latitude: "<<msg.latitude<<" Longitude: "<<msg.longitude<<" Altitude: "<<msg.altitude<<" WP"<<msg.index<<endl;			
	navicomout.publish(msg);

	while(ok()){
		
		if (mode == 1){
			if (mykbhit()){
				c = mygetch();
				switch (c){
					case '[' : Gas++; break;
					case ';' : Gas--; break;
					case 'd' : Roll=Roll+STEP; break;
					case 'a' : Roll=Roll-STEP; break;
					case 'w' : Pitch=Pitch+STEP; break;
					case 's' : Pitch=Pitch-STEP; break;
					case 'q' : Yaw=Yaw+STEP; break;
					case 'e' : Yaw=Yaw-STEP; break;
					case 'r' :
						{
							Gas = 126;
							Roll = 0;
							Pitch = 0;
							Yaw = 0;
						}break;
				}
				if (Gas<0) Gas=0;
			}
			cout<<"Key: "<<(char)i<<" Gas: "<<Gas<<" Roll: "<<Roll<<" Nickd: "<<Pitch<<" Yaw: "<<Yaw<<endl;
			mk_msgs::motorCommands msg;
			msg.throttle = Gas;
			msg.roll = Roll;
			msg.pitch = Pitch;
			msg.yaw = Yaw;		
			motorcomout.publish(msg);
		}else if (mode == 2){

		}
		Duration(0.1).sleep();
		spinOnce();
	}

	spin();
	return 0;
}
