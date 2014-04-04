/*
 * EarthCoords.h
 *
 *  Created on: Feb 12, 2011
 *      Author: joshua
 */

#ifndef EARTHCOORDS_H_
#define EARTHCOORDS_H_
#define _USE_MATH_DEFINES

#define MEAN_RADIUS_KM 6371
#define WORLD_POLAR_M 6356752.3142
#define WORLD_EQUATORIAL_M 6378137.0

#include <math.h>

namespace RSN{

/**
 * See :
 * http://www.movable-type.co.uk/scripts/latlong.html
 * for dicussion
 */
class EarthCoords {
public:
	EarthCoords();
	virtual ~EarthCoords();
	/**
	 * Returns the current home cords, or 0,0,0 if not set.
	 */
	void getHome(double&lat, double&lon, double&alt, double&roll ,double&pitch, double&yaw);
	/**
	 * Sets y = n/s with y+ --> north
	 * Sets x = e/w with x+ --> east
	 * Sets z = altitude with respect to the sea level 
	 * Note that, by necessity, all axis are Rhumb lines.
	 */
	void setHome(double lat, double lon, double alt, double roll, double pitch, double yaw);
	/**
	 * Returns x-y-z coordinates in the local world frame.
	 * x is e-w (longitude) distance in meters to base. y is n-s (latitude) distance in meters to base, z is difference in altitude
	 */
	void getLocalCords(double lat, double lon, double alt, double& xout, double& yout, double& zout);
	/**
	*/
	void getLocalAttit(double roll, double pitch, double yaw, double& rout,double& pout,double& yut);
	/**
	 * A spherical distance approximation is accurate ~ 1m. But it is fast.
	 */
	double getSphericalDistance(double lat1, double lon1, double lat2, double lon2);
	double degToRad(double);
private:
	bool homeset;
	double lathom;
	double lonhom;
	double althom;
	double rollhom;
	double pitchhom;
	double yawhom;
};
}

#endif /* EARTHCOORDS_H_ */
