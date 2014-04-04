/*
 * ClosedTrajectories.h
 *
 *  Created on: Jan 15, 2014
 *      Author: loki
 */

#ifndef CLOSEDTRAJECTORIES_H_
#define CLOSEDTRAJECTORIES_H_

namespace RSN {

class ClosedTrajectories {
protected:
	arma::vec3 center;
	arma::vec3 radius;
	arma::vec3 height;
public:
	ClosedTrajectories();
	virtual ~ClosedTrajectories();
	/**
	 * Here we specify where the center of the trajectory is located
	 */
	void setCenter(double x, double y, double z);
	/**
	 * This method sets the maximum distance from the center robot can travel
	 */
	void setRadius(double x, double y, double z);
	/**
	 * the altitude range is set by this method
	 * The idea is to give starting and the final altitude
	 */
	void setHeight(double low, double hight);
	/*
	 * Here the actual trajectory is generated...
	 */
	virtual void step();
};

} /* namespace RSN */

#endif /* CLOSEDTRAJECTORIES_H_ */
