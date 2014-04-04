/*
 * ClosedTrajectories.cpp
 *
 *  Created on: Jan 15, 2014
 *      Author: loki
 */

#include "ClosedTrajectories.h"

namespace RSN {

ClosedTrajectories::ClosedTrajectories() {
	//Initialize the variables
	center.fill(0.0);
	radius.fill(0.0);
	height.fill(0.0);
}

ClosedTrajectories::~ClosedTrajectories() {
}

void ClosedTrajectories::setCenter(double x, double y, double z){
	center(0) = x;
}


} /* namespace RSN */
