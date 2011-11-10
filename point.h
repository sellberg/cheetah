/*
 *  point.h
 *	-------
 *	Created by Jonas Sellberg on 2011-07-18.
 *
 *	Adapted from lecture notes of the course CS106B at Stanford University.
 *	-----------------------------------------------------------------------
 *	This interface provides a simple class that combines the
 *	coordinates of a point into a single object. The class Point
 *	exports a constructor for creating a new point of integers
 *	along with "getters" to extract the individual components.
 */

#ifndef _point_h
#define _point_h

class Point {
	
public:

	/*
	 * Constructor: Point(x, y)
	 * ------------------------
	 * Creates a new Point object with the indicated values.
	 */
	
	Point(int x, int y);
	
	/*
	 * Destructor: ~Point()
	 * --------------------
	 * Called when this Point is deleted or goes out of scope.
	 */
	
	~Point();
	
	/*
	 * Methods: getX(), getY()
	 * -----------------------
	 * These methods return the appropriate component.
	 */
	
	int getX();
	int getY();
		
private:
	
	/* Data required to implement a point */
	int x, y;
	
};

#endif
