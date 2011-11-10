/*
 *  peakdetect.cpp
 *	--------------
 *	Created by Jonas Sellberg on 2011-07-18.
 *	Copyright 2011 SLAC National Accelerator Laboratory. All rights reserved.
 *
 *  Object-oriented C++ version of the functionality available 
 *	in the MATLAB script at http://billauer.co.il/peakdet.html
 *	----------------------------------------------------------
 *	This file implements the PeakDetect class.
 */

#include "peakdetect.h"
#include "pointvector.h"
#include "point.h"

#include <iostream>

using namespace std;

/*
 * Implementation notes: PeakDetect constuctors
 * --------------------------------------------
 * The constructor must allocate the storage for the point
 * vectors of the maxima and minima and initialize the pointers
 * to the data values.
 */

PeakDetect::PeakDetect(uint16_t *Yarray, unsigned length) {
	x = NULL;
	y = Yarray;
	maxima = new PointVector();
	minima = new PointVector();
	this->length = length;
	index = 0;
}

PeakDetect::PeakDetect(int *Xarray, uint16_t *Yarray, unsigned length) {
	x = Xarray;
	y = Yarray;
	maxima = new PointVector();
	minima = new PointVector();
	this->length = length;
	index = 0;
}

/*
 * Implementation notes: ~PeakDetect
 * ---------------------------------
 * The destructor must deallocate the two point vectors
 */

PeakDetect::~PeakDetect() {
	delete maxima;
	delete minima;
}

/*
 * Implementation notes: clear
 * ---------------------------
 * This method clears all maxima and minima that has been found
 * and resets all instance variables.
 */

void PeakDetect::clear() {
	maxima->clear();
	minima->clear();
	index = 0;
}

/*
 * Implementation notes: findAll
 * -----------------------------
 * This method finds all maxima and minima that have a difference larger than delta.
 * It starts looking for a minima and adds an additional minima after the last maxima.
 */

void PeakDetect::findAll(float delta) {
	
	clear();
	
	Point max(-1, -1);
	Point min(-1, 65536);
	bool findMax = false;
	
	for (index = 0; index < length; index++) {
		// check if max/min should be updated
		if (y[index] > max.getY()) max = Point(index, y[index]);
		if (minima->size()) {
			if (y[index] < min.getY()) min = Point(index, y[index]);
		} else {
			if (y[index] <= min.getY()) min = Point(index, y[index]);
		}
		// check if max/min are larger than delta
		if (findMax) {
			if (y[index] < max.getY() - delta) {
				maxima->add(max.getX(), max.getY());
				min = Point(index, y[index]);
				findMax = false;
			}
		} else {
			if (y[index] > min.getY() + delta) {
				minima->add(min.getX(), min.getY());
				max = Point(index, y[index]);
				findMax = true;
			}
		}
	}
	// add last minimum
	minima->add(min.getX(), min.getY());
	
	if (x != NULL) {
		
		PointVector *maxtemp = new PointVector();
		PointVector *mintemp = new PointVector();
		
		for (int i = 0; i < maxima->size(); i++)
			maxtemp->add(x[(maxima->get(i))->getX()], (maxima->get(i))->getY());
		for (int i = 0; i < minima->size(); i++)
			mintemp->add(x[(minima->get(i))->getX()], (minima->get(i))->getY());
		
		delete maxima;
		delete minima;
		maxima = maxtemp;
		minima = mintemp;
		
	}
}

/*
 * Implementation notes: findNext
 * ------------------------------
 * This function 
 */

void PeakDetect::findNext(float delta) {

	Point max(-1, -1);
	Point min(-1, 65536);
	bool findMax = true;
	int maximaSize = 0;
	int minimaSize = 1;
	
	if (index) {
		max = Point(index, y[index]);
		min = max;
		maximaSize = maxima->size();
		minimaSize = minima->size();
	}
	
	while (maxima->size() - maximaSize < 1 || minima->size() - minimaSize < 1) {
		
		if (index >= length) break;
		// check if max/min should be updated
		if (y[index] > max.getY()) {
			if (x == NULL) max = Point(index, y[index]);
			else max = Point(x[index], y[index]);
		}
		if (y[index] < min.getY()) {
			if (x == NULL) min = Point(index, y[index]);
			else min = Point(x[index], y[index]);
		}
		// check if max/min are larger than delta
		if (findMax) {
			if (y[index] < max.getY() - delta) {
				maxima->add(max.getX(), max.getY());
				// find first minimum and add it
				if (maxima->size() == 1) {
					int jmax;
					if (x == NULL) jmax = max.getX();
					else {
						int k = 0;
						while (true) {
							if (x[k] == max.getX()) {
								jmax = k;
								break;
							}
							k++;
						}
					}
					min = Point(-1, 65536);
					for (int j = 0; j < jmax; j++) {
						if (y[j] <= min.getY()) {
							if (x == NULL) min = Point(j, y[j]);
							else min = Point(x[index], y[index]);
						}
					}
					minima->add(min.getX(), min.getY());
				}
				if (x == NULL) min = Point(index, y[index]);
				else min = Point(x[index], y[index]);
				findMax = false;
			}
		} else {
			if (y[index] > min.getY() + delta) {
				minima->add(min.getX(), min.getY());
				if (x == NULL) max = Point(index, y[index]);
				else max = Point(x[index], y[index]);
				findMax = true;
			}
		}
		if (++index == length) minima->add(min.getX(), min.getY());
	}
		
}
