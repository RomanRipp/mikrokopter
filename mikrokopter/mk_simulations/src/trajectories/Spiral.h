/*
 * Spiral.h
 *
 *  Created on: Jan 15, 2014
 *      Author: loki
 */

#ifndef SPIRAL_H_
#define SPIRAL_H_

class Spiral: public ClosedTrajectories {
public:
	Spiral();
	virtual ~Spiral();
	/**
	 * The single step on spiral trajectory
	 */
	void step();
};

#endif /* SPIRAL_H_ */

