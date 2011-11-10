/*
 *  pointvector.h
 *	------------
 *	Created by Jonas Sellberg on 2011-07-18.
 *
 *	Adapted from lecture notes of the course CS106B at Stanford University.
 *	-----------------------------------------------------------------------
 *	This interface defines the PointVector class, which implements
 *	a dynamic array of points made of integers.
 */

#ifndef _pointvector_h
#define _pointvector_h

#include "point.h"

/*
 * Class: PointVector
 * -----------------
 * This interface defines a class that models a vector of points.
 * The points are implemented by the Point class.
 * The fundamental point operations are add and get.
 */

class PointVector {
	
public:
	
	/*
	 * Constructor: PointVector
	 * Usage: PointVector pa;
	 * ----------------------
	 * Initializes a new empty vector that can contain points.
	 */
	
	PointVector();
	
	/*
	 * Destructor: ~PointVector
	 * Usage: (usually implicit)
	 * -------------------------
	 * Deallocates storage associated with this pa.  This method is
	 * called whenever a PointVector instance variable is deallocated.
	 */
	
	~PointVector();
	
	/*
	 * Method: size
	 * Usage: unsigned nElems = pa.size();
	 * --------------------------------
	 * Returns the number of points in this vector.
	 */
	
	unsigned size();
	
	/*
	 * Method: isEmpty
	 * Usage: if (pa.isEmpty()) . . .
	 * --------------------------------
	 * Returns true if this vector contains no points, and false otherwise.
	 */
	
	bool isEmpty();
	
	/*
	 * Method: clear
	 * Usage: pa.clear();
	 * --------------------
	 * This method removes all points from this vector.
	 */
	
	void clear();
	
	/*
	 * Method: add
	 * Usage: pa.add(x, y);
	 * ---------------------
	 * Adds the point defined by x and y onto this vector.
	 * ---------------------------------------------------
	 * Usage: pa.add(mypoint);
	 * ------------------------
	 * Adds the point defined by the pointer mypoint onto this vector.
	 */
	
	void add(int x, int y);
	void add(Point *point);
	
	/*
	 * Method: get
	 * Usage: Point *p = pa.get();
	 * -----------------------------
	 * Returns the pointer to the last point from this vector.
	 * -------------------------------------------------------
	 * Usage: Point *p = pa.get(index);
	 * ---------------------------------
	 * Returns the pointer to the point specified by index from this vector.
	 */
	
	Point *get();
	Point *get(unsigned index);	
	
private:
	
	/* Data required to implement a vector of points */
	
	Point **points;		/* Dynamic array of pointers to the points */
	unsigned capacity;  /* Allocated size of the array */
	unsigned count;     /* Number of points in the vector */
	
	/* Private method prototypes */
	
	void expandCapacity();
	
};

#endif
