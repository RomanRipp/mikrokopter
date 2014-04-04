/*
 * EarthCoords.cpp
 *
 *  Created on: Feb 12, 2011
 *      Author: joshua
 */

#include "EarthCoords/EarthCoords.h"
using namespace RSN;
double EarthCoords::degToRad(double deg){
	return (deg*M_PI/180.0);
}

void inline ellipticParameters(double a, double b, double phi,double&m, double&n){
		double eccentricity = acos(b/a);
		double n_prime = 1/(sqrt(1-pow(sin(phi), 2.0)*pow(sin(eccentricity), 2.0)));
		m = a*pow(cos(eccentricity), 2.0)*pow(n_prime, 3.0);
		n = a*n_prime;
}

EarthCoords::EarthCoords() {
	this->homeset=false;
}

EarthCoords::~EarthCoords() {

}
/**
 * Sets y = n/s with y+ --> north
 * Sets x = e/w with x+ --> east
 * Note that, by necessity, all axis are Rhumb lines.
 */
void EarthCoords::setHome(double lat, double lon){
	homeset = true;
	lathom=lat;
	lonhom=lon;
}
void EarthCoords::getHome(double&lat, double&lon){
	if (homeset){
		lat = lathom;
		lon=lonhom;
	} else {
		lat=0;
		lon=0;
	}
}
void EarthCoords::getLocalCords(double lat, double lon, double& xout, double& yout){
	double m,n;
	ellipticParameters(WORLD_EQUATORIAL_M,WORLD_POLAR_M,degToRad(lat),m,n);
	double diffLon = lon-lonhom;
	double diffLat = lat-lathom;
	double surfDistLong = M_PI/180.0*cos(degToRad(lat))*n;
	double surfDistLat = M_PI/180.0*m;
	xout = diffLon*surfDistLong;
	yout = diffLat*surfDistLat;
}
/**
 * A spherical distance approximation is accurate ~ 1m. But it is fast.
 */
double inline EarthCoords::getSphericalDistance(double lt1, double ln1, double lt2, double ln2){
	double R = MEAN_RADIUS_KM*1000.0; //we use mks round these parts.
	double dLat = degToRad(lt2-lt1);
	double dLon = degToRad(ln2-ln1);
	double a = sin(dLat/2) * sin(dLat/2) + cos(degToRad(lt1))
			* cos(degToRad(lt2))
			* sin(dLon/2) * sin(dLon/2);
	double c = 2 * atan2(sqrt(a), sqrt(1-a));
	double ret = R * c;
	return ret;
}
