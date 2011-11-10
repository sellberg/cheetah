/*
 *  peakdetect.h
 *	------------
 *	Created by Jonas Sellberg on 2011-07-18.
 *	Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
 *
 *  Object-oriented C++ version of the functionality available 
 *	in the MATLAB script at http://billauer.co.il/peakdet.html
 *	----------------------------------------------------------
 *	This interface defines the PeakDetect class, which implements
 *	peak detection for C-arrays.
 */

#ifndef _peakdetect_h
#define _peakdetect_h

#include "pointvector.h"

#include <stdint.h>

/*
 * Class: PeakDetect
 * -----------------
 * This interface defines a class that detects peaks in a C-array.

 * Characters are added and removed only from the top of the stack.
 * The fundamental stack operations are push (add to top) and pop
 * (remove from top).
 */

class PeakDetect {
	
public:
	
	/*
	 * Constructor: PeakDetect
	 * -----------------------
	 * Initializes a new PeakDetect from a C-array of unsigned 16-bit integers and its length.
	 * The C-array contains the Y-values for which the peak detection is performed.
	 * Assumes the corresponding X-value of each Y-value is determined by its index in the array.
	 */
	
	PeakDetect(uint16_t *Yarray, unsigned length);
	
	/*
	 * Constructor: PeakDetect
	 * -----------------------
	 * Initializes a new PeakDetect from two C-arrays of unsigned 16-bit integers and their length.
	 * The first C-array contains X-values.
	 * The second C-array contains Y-values.
	 * Prerequisites: Both C-arrays must have the same length!
	 */
	
	PeakDetect(int *Xarray, uint16_t *Yarray, unsigned length);
	
	/*
	 * Destructor: ~PeakDetect
	 * -----------------------
	 * Deallocates storage associated with this instance of PeakDetect.  This method is
	 * called whenever a PeakDetect instance variable is deallocated.
	 */
	
	~PeakDetect();
	
	/*
	 * Method: clear
	 * --------------------
	 * This method removes all maxima and minima from this PeakDetect and resets all instance variables.
	 */
	
	void clear();
	
	/*
	 * Method: findAll
	 * --------------------------------
	 * Finds all the peaks contained in the allocated Y-values whose heights are larger than delta.
	 */
	
	void findAll(float delta);
	
	/*
	 * Method: findNext
	 * --------------------------------
	 * Finds the next peak contained in the allocated Y-values whose height is larger than delta.
	 */
	
	void findNext(float delta);
	
	/* Public objects */
	
	PointVector *maxima;    // Pointer to a dynamic array of points holding the maxima
	PointVector *minima;    // Pointer to a dynamic array of points holding the minima
	
private:
	
	/* Data required to implement a PeakDetect */
	
	int *x;					// Pointer to array of X-values
	uint16_t *y;			// Pointer to array of Y-values
	unsigned length;		// Length of arrays of X- and Y-values
	unsigned index;			// Index for arrays of X- and Y-values
	
	/* Private method prototypes */
		
};

#endif
