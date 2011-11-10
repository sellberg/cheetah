/*
 *  point.cpp
 *	---------
 *	Created by Jonas Sellberg on 2011-07-18.
 *
 *	Adapted from lecture notes of the course CS106B at Stanford University.
 *	-----------------------------------------------------------------------
 *	This file implements the point.h interface.
 */

#include "point.h"

/*
 *	-------------------------------
 *	Constants
 *	-------------------------------
 */


/*
 *	-------------------------------
 *	Constructors and Destructor
 *	-------------------------------
 */

Point::Point(int x, int y) {
	this->x = x;
	this->y = y;
}

Point::~Point() {
	/* Empty */
}


/*
 *	-------------------------------
 *	Methods
 *	-------------------------------
 */

int Point::getX() {
	return x;
}

int Point::getY() {
	return y;
}
