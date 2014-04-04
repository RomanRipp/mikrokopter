#include <rendezvous/rendezvous.h>

#define LOOP_PERIOD_MS 50.0

using namespace std;

Rendezvous::Rendezvous(void) {
	setupROS();
}

Rendezvous::~Rendezvous(void) {
	nh.shutdown();
}

void Rendezvous::start(void) {
	
	ROS_INFO("Starting landing rendezvous");
	while(ros::ok()) {
		
		ros::Duration(LOOP_PERIOD_MS/1000.0).sleep();
		ros::spinOnce();
	}
}

void Rendezvous::setupROS(void) {
	
	//comunication service with mikrokopter and communiction client with husky
	interfacetouav_ser = nh.advertiseService(UAV_NAVI_DATA, &Rendezvous::getRendezvoiusPoint, this);
	interfacetougv_clt = nh.serviceClient<ugv_msgs::WaypointList>(UGV_NAVI_DATA, this);
	//Subscribe to husky current position and next waypoint
	currposugv_sub = nh.subscribe<geometry_msgs::Twist>(UGV_CURR_POSITION, 1000, &Rendezvous::getUgvCurrCoords, this);
	//Advertise rendevous point
	commonlat_pub = nh.advertise<std_msgs::Float32>(COMMON_LAT, 5, this);
	commonlon_pub = nh.advertise<std_msgs::Float32>(COMMON_LON, 5, this);
	//inintialize variables
	ugvCurrLat = 0;
	ugvCurrLon = 0;
	//ugvCurrElv = 0;
	uavCurrLat = 0;
	uavCurrLon = 0;
	commonLon = 0;
	commonLat = 0;
	//parameter initialization
	if (!nh.getParam("/rendezvous/parameter/ugvSpeed", Vg)){
		ROS_WARN("Cannot set UGV speed parameter, using default value 1 m/s");
		Vg = 2;
	} else {
		ROS_INFO("UGV speed: %f m/s", Vg);
	}
	if (!nh.getParam("/rendezvous/parameter/uavSpeed", Va)){
		ROS_WARN("Cannot set UAV speed parameter, using default value 2 m/s");
		Va = 4;
	} else {
		ROS_INFO("UAV speed: %f m/s", Va);
	}
}

void Rendezvous::getUgvCurrCoords(const geometry_msgs::TwistConstPtr &msg){
	ugvCurrLat = msg->linear.x;
	ugvCurrLon = msg->linear.y;
	//ugvCurrElv = msg->linear.z;
};

bool Rendezvous::getRendezvoiusPoint(mk_msgs::GetRendezvousPoint::Request& req,
									mk_msgs::GetRendezvousPoint::Response& res){

//	earthCoords.setHome(ugvCurrLat, ugvCurrLon);
	earthCoords.setHome(44.970969, -93.236648);
//	ROS_INFO("UGV current coordinates :%f, %f; UAV current coordinates: %f, %f", ugvCurrLat, ugvCurrLat, req.currLatitude, req.currLongitude);
	double x1,y1; //current position of ugv
	earthCoords.getLocalCords(44.970969, -93.236648,x1,y1);
	ROS_INFO("UGV current position: x: %f, y: %f",x1,y1);

	double x2,y2; // position of a ugv next waypoint
	//TODO get next waypoint coordinates
	earthCoords.getLocalCords(44.971086, -93.236884,x2,y2);
	ROS_INFO("UGV next waypoint position: x: %f, y: %f",x2,y2);

	double x3,y3; // position of a hexacopter
//	earthCoords.getLocalCords(req.currLatitude, req.currLongitude,x3,y3);
	earthCoords.getLocalCords(44.970786, -93.236707,x3,y3);
	ROS_INFO("UAV current position: x: %f, y: %f",x3,y3);

	double x4, y4;
	//TODO genegate landing point send and it back
	computePoint(x1,y1,x2,y2,x3,y3,x4,y4);
	earthCoords.getGlobalCords(x4,y4,commonLat,commonLon);
	ROS_INFO("Rendezvous point: x4: %f y4: %f lat: %f lon: %f", x4, y4, commonLat, commonLon);
	
/*
	ugv_msgs::WaypointList srv;
	ROS_INFO("Retrieving list of curent UGV waypoints");
//	srv.request.SET == 0;
//	srv.request.GET == 0;
	srv.request.command = 0; //command == GET - for returning WP list, have no idea why ?  
	if (!interfacetougv_clt.call(srv)){
		ROS_ERROR("Unable to retrieve current waypoint list of the UGV");
		return false;	
	}else{
		srv.response.cur_list;	
	}
	//TODO display waypoint list here

	ROS_INFO("Setting meeting point for UGV...");
//	srv.request.GET==0;
//	srv.request.SET==0;
	srv.request.command = 1;
//	srv.req.in_list.x = commonLat;
//	srv.req.in_list.y = commonLon;
//	srv.req.in_list.theta = ????;
	if (!interfacetougv_clt.call(srv)){
		ROS_ERROR("Unable to write rendezvous point to the UGV");
		return false;	
	}
*/
	ROS_INFO("Setting meeting point for UAV...");
	res.ack = 1;
	res.rendezLatitude = commonLat;
	res.rendezLongitude = commonLon;
	return true;
}

void Rendezvous::computePoint(double x1, double y1, double x2, double y2, double x3, double y3, double& x4, double& y4){
	ROS_INFO("Generating meeting point....");	
	if (Va==Vg){
    	double x4=(pow(x3,2)+pow(y3,2))/(2*(x3+(y2/x2)*y3));
    	double y4=(y2/x2)*x4;
	}else{
    	double a=(1+pow((y2/x2),2))*(pow(Vg,2)-pow(Va,2));
    	double b=(x3+(y2/x2)*y3)*(-2)*pow(Vg,2);
    	double c=(pow(x3,2)+pow(y3,2))*pow(Vg,2);
    	double D=pow(b,2)-4*a*c;
    	double X4[2] = {(-b+sqrt(D))/(2*a), (-b-sqrt(D))/(2*a)};
    	double Y4[2] = {(y2/x2)*X4[0], (y2/x2)*X4[1]};
        
    	int eq1= (y2/x2)*X4[0]-Y4[0] + sqrt(pow(X4[0],2)+pow(Y4[0],2))/Vg-sqrt(pow((X4[0]-x3),2)+pow((Y4[0]-y3),2))/Va;
    	int eq2= (y2/x2)*X4[1]-Y4[1] + sqrt(pow(X4[1],2)+pow(Y4[1],2))/Vg-sqrt(pow((X4[1]-x3),2)+pow((Y4[1]-y3),2))/Va;
    	if (eq1 > 0 || eq2 > 0)
        	ROS_ERROR("Invalid roots!");
    	if (((X4[0]>0 && X4[0]<x2) || (X4[0]<0 && X4[0]>x2)) || ((Y4[0]>0 && Y4[0]<y2) || (Y4[0]<0 && Y4[0]>y2))){
        	x4 = X4[0];
        	y4 = Y4[0];
    	}else if (((X4[1]>0 && X4[1]<x2) || (X4[1]<0 && X4[1]>x2)) || ((Y4[1]>0 && Y4[1]<y2) || (Y4[1]<0 && Y4[1]>y2))){
            	x4 = X4[1];
            	y4 = Y4[1];
        	}else{
            	ROS_INFO("Lets meet at UGV current position");
            	x4 = 0;
            	y4 = 0;
			}
    }
}

int main(int argc, char** argv) {

	ros::init(argc,argv,"mk_rendezvous");
	Rendezvous r;

	r.start(); // blocking

	return 0;
}
