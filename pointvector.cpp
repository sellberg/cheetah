/*
 *  pointvector.cpp
 *	---------------
 *	Created by Jonas Sellberg on 2011-07-18.
 *
 *	Adapted from lecture notes of the course CS106B at Stanford University.
 *	-----------------------------------------------------------------------
 *	This file implements the PointVector class.
 */

#include "pointvector.h"
#include "point.h"

#include <iostream>

/*
 * Constants: INITIAL_CAPACITY
 * --------------------------------------
 * This constant defines the initial allocated size of the points
 * array.  If the vector grows beyond its capacity, the implementation
 * doubles the allocated size.
 */

const int INITIAL_CAPACITY = 100;

/*
 * Implementation notes: PointVector constuctor
 * ------------------------------------------
 * The constructor must allocate the array storage for the vector
 * points and initialize the fields of the object.
 */

PointVector::PointVector() {
	capacity = INITIAL_CAPACITY;
	points = new Point *[capacity];
	for (int i = 0; i < capacity; i++) {
		points[i] = NULL;
	}	
	count = 0;
}

/*
 * Implementation notes: ~PointVector
 * --------------------------------
 * The destructor must deallocate every point (which it can do by
 * calling clear) and then free the dynamic array.
 */

PointVector::~PointVector() {
	clear();
	delete[] points;
}

/*
 * Implementation notes: size, isEmpty
 * -----------------------------------
 * These implementations should be self-explanatory.
 */

unsigned PointVector::size() {
	return count;
}

bool PointVector::isEmpty() {
	return count;
}

/*
 * Implementation notes: clear
 * ---------------------------
 * This method deallocates every point that has been allocated
 * in the dynamic array.
 */

void PointVector::clear() {
	for (int i = 0; i < count; i++) {
		delete points[i];
		points[i] = NULL;
	}
	count = 0;
}

/*
 * Implementation notes: add
 * --------------------------
 * This function must first check to see whether there is
 * enough room for the point and expand the array storage
 * if necessary. If safe, it adds a new Point to the heap
 * or copies the address of an already existing Point.
 */

void PointVector::add(int x, int y) {
	if (count == capacity) expandCapacity();
	points[count++] = new Point(x, y);
}

void PointVector::add(Point *point) {
	if (count == capacity) expandCapacity();
	points[count++] = point;
}

/*
 * Implementation notes: get
 * -------------------------------
 * These functions must check for an empty vector and report
 * NULL if there is no point at the vector index specified.
 */

Point *PointVector::get() {
	if (isEmpty()) return NULL;
	return points[count - 1];
}

Point *PointVector::get(unsigned index) {
	if (index < 0 || index > count - 1) return NULL;
	return points[index];
}

/*
 * Implementation notes: expandCapacity
 * ------------------------------------
 * This private method doubles the capacity of the points array
 * whenever it runs out of space.  To do so, it must allocate a new
 * array, copy all the points from the old array to the new one,
 * free the old storage, and finally replace the points pointer
 * with the new array.
 */

void PointVector::expandCapacity() {
	capacity *= 2;
	Point **array = new Point *[capacity];
	for (int i = 0; i < capacity; i++) {
		if (i < count) array[i] = points[i];
		else array[i] = NULL;
	}	
	delete[] points;
	points = array;
}
