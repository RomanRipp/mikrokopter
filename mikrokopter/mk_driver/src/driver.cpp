#include "driver.h"
#include <boost/thread.hpp>

using namespace std;
using namespace ros;

Publisher sensordata_pub;		        
Publisher statedata_pub;

speed_t getBaud(int speed){
	switch(speed){
		case 9600:
         		return B9600;
      		case 19200:
         		return B19200;
       		case 38400:
				return B38400;
	       	case 57600:
         		return B57600;
       		case 115200:
         		return B115200;
       		case 230400:
         		return B230400;
      		default:                 // invalid bitrate
         		return B0; 	
	} 
}

int openSerial(char *portname, int speed){
	ROS_INFO("Opening serial port for communication");

	struct termios options;
	tcgetattr(fd, &options);
	speed_t serialBaud = getBaud(speed);
	cfsetispeed(&options, serialBaud);
	cfsetospeed(&options, serialBaud);
	tcsetattr(fd, TCSANOW, &options);

	fd = open (portname, O_RDWR | O_NOCTTY | O_NDELAY);  
  	if (fd < 0)
  	{
		ROS_INFO("Error opening serial port");
		return -1;
  	}else
  	{
		fcntl(fd, F_SETFL, 0);
		ROS_INFO("Serial port opened");
  	}

	options.c_cflag &= ~PARENB;    // set no parity, stop bits, data bits
	options.c_cflag &= ~CSTOPB;
	options.c_cflag &= ~CSIZE;
	options.c_cflag |= CS8;
	tcsetattr(fd, TCSANOW, &options);
	return 1;
}

int receiveData(){
	unsigned char buf[RXD_BUFFER_LENGTH];
	int ReceivedBytes = 0;
	int ptr_rxd_buffer = 0;
	int nLength = 0;
	unsigned char c;
	unsigned char crc1, crc2;
	unsigned int crc;
	int l=0;

	// Initialize file descriptor sets
	fd_set read_fds, write_fds, except_fds;
	FD_ZERO(&read_fds);
	FD_ZERO(&write_fds);
	FD_ZERO(&except_fds);
	FD_SET(fd, &read_fds);

	// Set timeout to 1.0 seconds
	struct timeval timeout;
	timeout.tv_sec = LINK_TIMEOUT_S;
	timeout.tv_usec = LINK_TIMEOUT_MICROS;

	// Wait for input to become ready or until the time out; the first parameter is
	// 1 more than the largest file descriptor in any of the sets
	if (select(fd + 1, &read_fds, &write_fds, &except_fds, &timeout) == 1)
	{
    	// fd is ready for reading	
		nLength = read(fd, buf, sizeof buf);
	} else {
    	// timeout or error
		ROS_ERROR("Serial link timeout");
		unsigned char interval = NAVI_SEND_INTERVAL; //1 byte sending interval ( in 10ms steps ) 
		encode('o',0,&interval, 1);
		return -1;
	}
	while(buf[l]!='#' && l<=RXD_BUFFER_LENGTH)
		{
			l++;
		}
	for (int i = l; i<nLength; i++){
		c = buf[i];
		if ((ptr_rxd_buffer == 0) && (c == '#')){
			g_rxd_buffer[ptr_rxd_buffer++]=c;
			crc=c;		
		}
		else if (ptr_rxd_buffer < RXD_BUFFER_LENGTH){
			if ((int) c != 10){		                                 					
				g_rxd_buffer[ptr_rxd_buffer++] = c; // copy byte to rxd buffer
				crc += c; 
                         }
                         else // termination character was received
                         {
                crc -= g_rxd_buffer[ptr_rxd_buffer-2];
                crc -= g_rxd_buffer[ptr_rxd_buffer-1];
                // calculate checksum from transmitted data
                crc %= 4096;
                crc1 = '=' + crc / 64;
				crc2 = '=' + crc % 64;
				if((crc1 == g_rxd_buffer[ptr_rxd_buffer-2]) && (crc2 == g_rxd_buffer[ptr_rxd_buffer-1]))
				{   // checksum valid
					g_rxd_buffer[ptr_rxd_buffer] = '\r';
					ReceivedBytes = ptr_rxd_buffer + 1;
				}
				else
				{       // checksum invalid
					ptr_rxd_buffer = 0; 
					return -1;
				}
				ptr_rxd_buffer = 0; 
				}
			}
			else 
			{
				ptr_rxd_buffer = 0;
			}
	}
return ReceivedBytes;
} 

void decode(int receivedBytes) {
	unsigned char a,b,c,d;
    unsigned char x,y,z;
    int ptrIn = 3;
    int ptrOut = 3;
	int len = receivedBytes;
	int i=0;
 
         while(len)
         {	
                 a = g_rxd_buffer[ptrIn++] - '=';
                 b = g_rxd_buffer[ptrIn++] - '=';
                 c = g_rxd_buffer[ptrIn++] - '=';
                 d = g_rxd_buffer[ptrIn++] - '=';

                 x = (a << 2) | (b >> 4);
                 y = ((b & 0x0f) << 4) | (c >> 2);
                 z = ((c & 0x03) << 6) | d;

                 if(len--) {
				i=ptrOut++;
				g_rxd_buffer[i] = x; 
		} else break;
                 if(len--) {
				i=ptrOut++;
				g_rxd_buffer[i] = y; 
		} else break;
                 if(len--) {
				i=ptrOut++;
				g_rxd_buffer[i] = z; 
		} else break;
         }
	 g_pRxData =(unsigned char*) &g_rxd_buffer[3];
         g_RxDataLen = ptrOut - 3;
}

void publish(){
	switch(g_rxd_buffer[2]){
		case 'O':{
			//Navigation data
			NaviData_t NaviData;	
			memcpy((unsigned char *)&NaviData,(unsigned char *) g_pRxData,sizeof(NaviData));

			float lat, lon, alt, roll, pitch, yaw, heading, compassheading;
			lat = NaviData.CurrentPosition.Latitude / 10000000.0f;
			lon = NaviData.CurrentPosition.Longitude / 10000000.0f;
			alt = NaviData.CurrentPosition.Altitude / 1000.0f;
			roll = MIKO_TO_RAD(NaviData.AngleRoll);
			pitch = MIKO_TO_RAD(NaviData.AngleNick);
			heading = DEG_TO_RAD(NaviData.Heading);
			compassheading = DEG_TO_RAD(NaviData.CompassHeading);

			if (compassheading <= PI) {
				yaw = (-1)*DEG_TO_RAD(NaviData.CompassHeading);
			} else {
				yaw = 2*PI - DEG_TO_RAD(NaviData.CompassHeading);
			}
			if (NaviData.SatsInUse >= 4 && NaviData.FlyingTime == 0) {
				earthCoords.setHome(lat, lon, alt, roll, pitch, PI);
			}else if (NaviData.SatsInUse >= 4 && NaviData.FlyingTime > 0){
				double x,y,z,wroll,wpitch,wyaw;
				mk_msgs::stateData data;
				earthCoords.getLocalCords(lat, lon, alt, x, y, z);
				earthCoords.getLocalAttit(roll, pitch, yaw, wroll, wpitch, wyaw);
				data.x=x;
				data.y=y;
				data.z=z;
				data.roll=wroll;
				data.pitch=wpitch;
				data.yaw=yaw;
				statedata_pub.publish(data);
			}

			mk_msgs::sensorData buffer;
			buffer.Longitude=lon; 
			buffer.Latitude=lat;
			buffer.gpsAltitude=alt;
			buffer.NumberOfWPs=NaviData.WaypointNumber; 		//number of waypoints			
			buffer.CurrentWP=NaviData.WaypointIndex; 			//Index of a current Waypoint
			buffer.SatsInUse=NaviData.SatsInUse;				// number of satellites used for position solution
			buffer.Altitude=NaviData.Altimeter/10.0f;			//Barometric altitude in meters	
			buffer.FlyingTime=NaviData.FlyingTime;				// in seconds
			buffer.Battery=NaviData.UBat/10.0f;					// Battery Voltage in Volts		
			buffer.GroundSpeed=NaviData.GroundSpeed/100.0f;		// speed over ground in m/s (2D)			
			buffer.Heading=heading;								// current flight direction in RAD as angle to north
			buffer.CompassHeading=yaw;						    // current compass value in RAD
			buffer.Nick=pitch;									// current Nick angle in radians
			buffer.Roll=roll;									// current Nick angle in radians
			buffer.RCQuality=NaviData.RC_Quality;				// RC_Quality
			buffer.ErrorCode=NaviData.Errorcode;				// 0->OK
			buffer.zSpeed=NaviData.TopSpeed/100.0f;            	// velocity in vertical direction in m/s
			buffer.TargetHoldTime=NaviData.TargetHoldTime; 		// time in s to stay at the given target, counts down to 0 if target has been reached
	        buffer.Gas=NaviData.Gas;	         				// for future use
	        buffer.Current=NaviData.Current/10.0f;				// actual current in A steps
	        buffer.UsedCapacity=NaviData.UsedCapacity;			// used capacity in mAh
//			buffer.NCMode = NaviData.NCFlags;					//CH PH FREE and other usefull stuff
//			buffer.FCMode1 = NaviData.FCStatusFlags;			//low battery and other usefull stuff
//			buffer.FCMode2 = NaviData.FCStatusFlags2;			//Carefree alt hold and other usefull stuff
			
//			ROS_Flags("INFO: %d, TargetReached: %d",NaviData.NCFlags,NaviData.NCFlags & 0x20);
		
			buffer.Free = (NaviData.NCFlags & NC_FLAG_FREE);
			buffer.PH = (NaviData.NCFlags & NC_FLAG_PH)/NC_FLAG_PH;
			buffer.CH = (NaviData.NCFlags & NC_FLAG_CH)/NC_FLAG_CH;
			buffer.RangeLimit = (NaviData.NCFlags & NC_FLAG_RANGE_LIMIT)/NC_FLAG_RANGE_LIMIT;
			buffer.NoSerialLink = (NaviData.NCFlags & NC_FLAG_NOSERIALLINK)/NC_FLAG_NOSERIALLINK;
			buffer.TargetReached = (NaviData.NCFlags & NC_FLAG_TARGET_REACHED)/NC_FLAG_TARGET_REACHED;
			buffer.Manual = (NaviData.NCFlags & NC_FLAG_MANUAL)/NC_FLAG_MANUAL;
			buffer.GPSOK = (NaviData.NCFlags & NC_FLAG_GPS_OK)/NC_FLAG_GPS_OK;

			buffer.MotorsOn = (NaviData.FCStatusFlags & FC_STATUS_MOTOR_RUN)/FC_STATUS_MOTOR_RUN;
			buffer.Flying = (NaviData.FCStatusFlags & FC_STATUS_FLY)/FC_STATUS_FLY;
			buffer.LowBat = (NaviData.FCStatusFlags & FC_STATUS_LOWBAT)/FC_STATUS_LOWBAT;

			buffer.CareFree = (NaviData.FCStatusFlags2 & FC_STATUS2_CAREFREE_ACTIVE)/FC_STATUS2_CAREFREE_ACTIVE;
			buffer.AltHld = (NaviData.FCStatusFlags2 & FC_STATUS2_ALTITUDE_CONTROL_ACTIVE)/FC_STATUS2_ALTITUDE_CONTROL_ACTIVE;
			buffer.Failsafe = (NaviData.FCStatusFlags2 & FC_STATUS2_FAILSAFE_ACTIVE)/FC_STATUS2_FAILSAFE_ACTIVE;

			sensordata_pub.publish(buffer);	
		}			
		break;
		case 'P':{ 
			ROS_INFO("Serial Potis:");
			cout<<"Poti1: "<<(int)g_rxd_buffer[3]<<endl;
		}		
		break;

	}
}

int rescale(int b1, int b2){
	int val = 0;
	if(b2<122) val = b1+b2*255;
	else val = (b1-255)+(b2-255)*255;
	return val;
}

void motorCallback(const mk_msgs::motorCommands msg){	
	
	str_ExternControl ExternControl;
	ExternControl.Throttle = GAS_COEFF * msg.throttle;
	ExternControl.Roll = RAD_TO_MIKO(msg.roll);
    ExternControl.Pitch = RAD_TO_MIKO(msg.pitch);
    ExternControl.Yaw = RAD_TO_MIKO(msg.yaw);
	ExternControl.Config =1;
	ExternControl.Frame ='f';
	encode('b', 1, (unsigned char*) &ExternControl, sizeof(ExternControl));
}

bool setwpCallback(
			mk_msgs::SetGPSWaypoint::Request  &req,
         	mk_msgs::SetGPSWaypoint::Response &res){

	Point_t pointSent;
	pointSent.Position.Longitude = req.longitude * 10000000;	
	pointSent.Position.Latitude = req.latitude * 10000000; 
	pointSent.Position.Altitude = req.altitude * 10;
	pointSent.Position.Status = req.status;
	pointSent.HoldTime = req.holdTime;
	pointSent.Index = req.index;
	pointSent.AltitudeRate = req.altitudeRate*10;
	pointSent.Speed = req.speed*10;

	pointSent.ToleranceRadius = GPS_TOLERANCE_RADIUS;
	pointSent.Heading = 1;
	pointSent.Event_Flag = 0;
	pointSent.Type = 0;
	pointSent.WP_EventChannelValue = 0;
	pointSent.CamAngle = 0;
	pointSent.Name[0] = 'P';
	pointSent.Name[1] = '0'+req.index;
	pointSent.Name[2] = 0;
	pointSent.Name[3] = 0;
	pointSent.reserve[0] = 0;
	pointSent.reserve[1] = 0;
	ROS_INFO("Writing waypoint %d", pointSent.Index);
	encode('w', 2, (unsigned char*) &pointSent, sizeof(pointSent));

	int St = 0;
	int timeout = 0;
	while (res.ack != 1) {
		//Receive from serial
		St = receiveData();
		timeout++;
		if (St>0 && g_rxd_buffer[2] == 'W'){
			decode(St);
			//Waypoint confirmation 
			if (pointSent.Index == (int)g_rxd_buffer[3]) {
				ROS_INFO("Waypoint checked");
				res.ack = 1;
			}
		}else if (timeout*DRIVER_LOOP_MS >= RESEND_INTERVAL_MS){
			ROS_WARN("Waypoint write timeout occured, resending...");
			timeout = 0;
			ROS_INFO("Writing waypoint %d", pointSent.Index);
			encode('w', 2, (unsigned char*) &pointSent, sizeof(pointSent));
		}
		
	Duration(DRIVER_LOOP_MS/1000.0).sleep();
	}
	return true;
}

bool listwpCallback(
			mk_msgs::listwp::Request  &req,
         	mk_msgs::listwp::Response &res){
	int numOfWP = 1;
	int St = 0;
	int timeout = 0;
	ROS_INFO("Listing all the waypoints: "); 
	for (int i = 1; i <= numOfWP; i++){
		unsigned char buffer = i;
		encode('x', 2, &buffer, 1);
		bool printed = false;
		while (!printed) {
			//Receive from serial
			St = receiveData();
			timeout++;
			if (St>0 && g_rxd_buffer[2] == 'X'){
				decode(St); 
				numOfWP = (int)g_rxd_buffer[3];
				ROS_INFO("Number of waypoints: %d, Index: %d", numOfWP, i);

				g_pRxData =(unsigned char*) &g_pRxData[2];
				Point_t pointReceived;
				memcpy((unsigned char *)&pointReceived, (unsigned char *)g_pRxData, sizeof(pointReceived));
				ROS_INFO("Index      %d", (int)pointReceived.Index);
				ROS_INFO("Latitude:  %f", (int)pointReceived.Position.Latitude / 10000000.0f);
				ROS_INFO("Longitude: %f", (int)pointReceived.Position.Longitude / 10000000.0f);
				ROS_INFO("Altitude:  %f", (int)pointReceived.Position.Altitude / 10.0f);
				ROS_INFO("Status:    %d", (int)pointReceived.Position.Status);
				ROS_INFO("AltRate:   %f", (int)pointReceived.AltitudeRate / 10.0f);
				ROS_INFO("Heading:   %d", (int)pointReceived.Heading);
				ROS_INFO("HoldTime:  %d", (int)pointReceived.HoldTime);
				ROS_INFO("Speed:     %f", (int)pointReceived.Speed / 10.0f);
				ROS_INFO("Radius:    %d", (int)pointReceived.ToleranceRadius);
				printed = true;
			}else if (timeout*DRIVER_LOOP_MS == RESEND_INTERVAL_MS){
				ROS_WARN("Waypoint list timeout occured, resending...");
				timeout = 0;
				encode('x', 2, &buffer, 1);
			}
		Duration(DRIVER_LOOP_MS/1000.0).sleep();
		}
	}
	res.ack = 1;
	return true;
}

bool setptiCallback(
			mk_msgs::SetSerialPoti::Request  &req,
         	mk_msgs::SetSerialPoti::Response &res){
	
	unsigned char* potis = new unsigned char[12];
	potis[0] = req.poti1;
	potis[1] = req.poti2;
	potis[2] = req.poti3;
	potis[3] = req.poti4;
	for(int i = 4; i<12; i++) potis[i] = -125;
	ROS_INFO("Switching mikrokopter flight mode...");

	encode('y', 1, potis, 12);
	unsigned char interval = NAVI_SEND_INTERVAL; //1 byte sending interval ( in 10ms steps ) 
	encode('o',0,&interval, 1);

	int St = 0;
	int timeout = 0;
	ROS_INFO("Checking if mode switched correcty...");
	while (res.ack != 1) {
		//Receive from serial
		St = receiveData();
		timeout++;
		if (St>0 && g_rxd_buffer[2] == 'O'){
			decode(St);
			if (ModeSwitched(req)){
				//Everything is fine
				res.ack = 1;
			}else{
				
			}
		}else if (timeout*DRIVER_LOOP_MS >= RESEND_INTERVAL_MS){
				ROS_WARN("Set poti timeout occured, resending...");
				timeout = 0;
				encode('y', 1, potis, 12);
				encode('o',0,&interval, 1);
		}
	Duration(DRIVER_LOOP_MS/1000.0).sleep();
	}
	return true; 
}

bool ModeSwitched(mk_msgs::SetSerialPoti::Request potis){

	bool switched = false;
	NaviData_t NaviData;	
	memcpy((unsigned char *)&NaviData,(unsigned char *) g_pRxData,sizeof(NaviData));
	int carefree = (int)NaviData.FCStatusFlags2 & FC_STATUS2_CAREFREE_ACTIVE;
	int gps = (int)NaviData.NCFlags;
	int althold = (int)NaviData.FCStatusFlags2 & FC_STATUS2_ALTITUDE_CONTROL_ACTIVE;
	if ((potis.poti1 == POTI_OFF && carefree == 0) 
	|| (potis.poti1 == POTI_ON && carefree == FC_STATUS2_CAREFREE_ACTIVE)){
		ROS_INFO("Carefree ok");
		if ((potis.poti2 == POTI_OFF && ((gps & NC_FLAG_FREE) == NC_FLAG_FREE)) 
		|| (potis.poti2 == POTI_NEUTRAL && ((gps & NC_FLAG_PH) == NC_FLAG_PH)) 
		|| (potis.poti2 == POTI_ON && ((gps & NC_FLAG_CH) == NC_FLAG_CH))){
			ROS_INFO("GPS ok");
			if ((potis.poti3 == POTI_OFF && althold == 0) 
			|| (potis.poti3 == POTI_ON && althold == FC_STATUS2_ALTITUDE_CONTROL_ACTIVE)){
				ROS_INFO("Alt hold ok");
				ROS_INFO("Mode changed to:");
				ROS_INFO("Carefree: %d",carefree);
				ROS_INFO("GPS:FREE: %d, PH: %d, CH: %d", gps & NC_FLAG_FREE, (gps & NC_FLAG_PH) / NC_FLAG_PH, (gps & NC_FLAG_CH) / NC_FLAG_CH);
				ROS_INFO("Alt hold: %d",althold/FC_STATUS2_ALTITUDE_CONTROL_ACTIVE);
				ROS_INFO("Ext ctrl: %d",(potis.poti4 + 125)/250);
				switched = true;
			} else {switched = false;}
		} else {switched = false;}
	} else {switched = false;}
	return switched;
}

void encode(unsigned char cmd,int addr, unsigned char* data, int dataLen)
{
	unsigned char txd_buffer[dataLen+1];
	unsigned int pt = 0;
	unsigned char a,b,c;
	unsigned char ptr = 0;
	txd_buffer[pt++] = '#';               // Start-Byte
	txd_buffer[pt++] = 'a' + addr;        // Adress
     	txd_buffer[pt++] = cmd;               // Command
     	while(dataLen)
     	{
      	if(dataLen) { a = data[ptr++]; dataLen--;} else a = 0;
      	if(dataLen) { b = data[ptr++]; dataLen--;} else b = 0;
      	if(dataLen) { c = data[ptr++]; dataLen--;} else c = 0;
  	txd_buffer[pt++] = '=' + (a >> 2);
        txd_buffer[pt++] = '=' + (((a & 0x03) << 4) | ((b & 0xf0) >> 4));
        txd_buffer[pt++] = '=' + (((b & 0x0f) << 2) | ((c & 0xc0) >> 6));
        txd_buffer[pt++] = '=' + ( c & 0x3f);
       	}	

	unsigned int tmpCRC = 0;
	int i = 0;
	for(i = 0; i <  pt; i++)
	{
		tmpCRC += txd_buffer[i];
	}
	tmpCRC %= 4096;
	txd_buffer[i++] = '=' + tmpCRC / 64;
	txd_buffer[i++] = '=' + tmpCRC % 64;
	txd_buffer[i++] = '\r';

	write(fd, txd_buffer, i);
}

int validSpeed(char* in){
	while (true){
		int speed = atoi(in);
		if (speed == 9600 ||
 			speed == 19200 || 
			speed == 38400 || 
			speed == 57600 || 
			speed == 115200 || 
			speed == 230400){		
			return speed;
		}else{
			cout<<"Baud speed value is incorrect, try again"<<endl;
			cin>>in;
		}
	}
return 0;
}

int main(int argc, char** argv){
  	ROS_INFO("Hexacopter wrapper started...");
	int baudSp = validSpeed(argv[2]);
	int St = 0; 
	unsigned char buffer=0; 

  	init(argc,argv,"wrapper");
  	NodeHandle nh("/mikrokopter/");

	//control structure initialization
	Subscriber motor_sub = nh.subscribe<mk_msgs::motorCommands>("commands/motor", 1000, motorCallback);
	ServiceServer setwp_srv = nh.advertiseService("commands/setwp", setwpCallback);
	ServiceServer listwp_srv = nh.advertiseService("commands/listwp", listwpCallback);
	ServiceServer setpti_srv = nh.advertiseService("commands/setpoti", setptiCallback);

	//Topic initialization v
	sensordata_pub = nh.advertise<mk_msgs::sensorData>("sensor/data", 1000);		        
	statedata_pub = nh.advertise<mk_msgs::stateData>("state/data", 1000);

	if (openSerial(argv[1], baudSp) < 0) return -1;

	//Navi Error request
	encode('e',2,0, 0);

	int count = 70;
  	Rate Loop_rate(1000.0/DRIVER_LOOP_MS);
  	while(ok()){
		count++;
		St = receiveData(); //receiving data from FC
		if (St>0){
			decode(St);
			publish();
		}
		//OSD Data request
		if (count >= 70){
			buffer = NAVI_SEND_INTERVAL; //1 byte sending interval ( in 10ms steps ) 
			encode('o',0,&buffer, 1);
			count = 0;
		}

		spinOnce();
		Loop_rate.sleep();
  	}
 	spin();
  	return 0;
}
