#ifndef RENDEZVOUS_H
#define RENDEZVOUS_H

#include <ros/ros.h>
#include <rendezvous/topicNames.h>
#include <mk_msgs/GetRendezvousPoint.h>
#include <ugv_msgs/WaypointList.h>
#include <geometry_msgs/Twist.h>
#include <std_msgs/Float32.h>
#include <EarthCoords/EarthCoords.h>

class Rendezvous {

	public:
		/**
		 * Constructor
		 */
		Rendezvous();

		/**
		 * Destructor.
		 */
		~Rendezvous();

		/**
		 * Starts the rendezvous
		 */
		void start(void);
	private:
	ros::NodeHandle nh;
	ros::ServiceServer interfacetouav_ser;
	ros::ServiceClient interfacetougv_clt;
	ros::Publisher commonlat_pub, commonlon_pub;
	ros::Subscriber currposugv_sub;
	RSN::EarthCoords earthCoords;
	double ugvCurrLat, ugvCurrLon, uavCurrLat, uavCurrLon, commonLat, commonLon, Vg, Va;//,ugvCurrElv;

	void setupROS(void);
	bool getRendezvoiusPoint(mk_msgs::GetRendezvousPoint::Request& req,	
							mk_msgs::GetRendezvousPoint::Response& res);
	void getUgvCurrCoords(const geometry_msgs::TwistConstPtr &msg);
	void computePoint(double x1, double y1, double x2, double y2, double x3, double y3, double& x4, double& y4);
};
#endif /* RENDEZVOUS_H */
