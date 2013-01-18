/*
 *  attenuation.cpp
 *
 *	Written by Jonas Sellberg on 2011-04-20 as a module to
 *
 *  cheetah
 *
 *  Created by Anton Barty on 6/2/11.
 *  Copyright 2011 CFEL. All rights reserved.
 *
 *You can modify this software under the terms of the GNU General Public License 
 *as published by the Free Software Foundation, either version 3 of the License, or
 *(at your option) any later version.
 *
 *This program is distributed in the hope that it will be useful,
 *but WITHOUT ANY WARRANTY; without even the implied warranty of
 *MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *GNU General Public License for more details.
 *
 *You should have received a copy of the GNU General Public License
 *along with this program.  If not, see <http://www.gnu.org/licenses/>. 
 */

#include "myana/myana.hh"
#include "myana/main.hh"
#include "release/pdsdata/cspad/ConfigV1.hh"
#include "release/pdsdata/cspad/ConfigV2.hh"
#include "release/pdsdata/cspad/ConfigV3.hh"
#include "release/pdsdata/cspad/ConfigV4.hh"
#include "release/pdsdata/cspad/ElementHeader.hh"
#include "release/pdsdata/cspad/ElementIterator.hh"

#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <iostream>
#include <string>
using std::cout;
using std::endl;
using std::string;

#include <fstream>
#include <sstream>

#include "setup.h"
#include "worker.h"
#include "attenuation.h"


/*
 *	Functions called by setup.cpp and cheetah.cpp to calculate attenuation
 */

// Function to calculate the Si thickness inside the beam
int getSiThickness(unsigned & totalThickness, unsigned nFilters, unsigned filterThicknesses[])
{
	int fail = 0;
	int totalFail = 0;
	float filterPositions[nFilters];
	int filterNumber;
	string pvName;
	totalThickness = 0;
	for (int i=0; i<nFilters; i++) {
		filterNumber = 2+i;
		if (filterNumber < 10) pvName = "XRT:DIA:MMS:0" + intToString(filterNumber) + ".RBV";
		else pvName = "XRT:DIA:MMS:" + intToString(filterNumber) + ".RBV";
		fail = getPvFloat(pvName.c_str(), filterPositions[i]);
		totalFail += fail;
		if (fail == 0) {
			// cout << "filterPosition: " << filterPositions[i] << endl;
			if (filterPositions[i] > 20 && filterPositions[i] < 30) totalThickness += filterThicknesses[i];
			else if (filterPositions[i] < -2.5 || filterPositions[i] > 2.5) {
				cout << "Silicon Filter " << filterNumber << " with thickness " << filterThicknesses[i] << " um is out of position: " << filterPositions[i] << " mm" << endl;
				totalFail++;
			}
		} else {
			cout << "Could not retrieve Silicon Filter " << (filterNumber-1) << " with thickness " << filterThicknesses[i] << " um" << endl;
			// cout << "filterPosition: " << filterPositions[i] << endl;
		}
	}
	return totalFail;
}


// Function to calcualte the atteunation of Si from a specific thickness
double getAttenuation(unsigned totalThickness, unsigned nThicknesses, unsigned possibleThicknesses[], double possibleAttenuations[])
{
	for (int i=0; i<nThicknesses; i++) {
		if (possibleThicknesses[i] == totalThickness) return 1/possibleAttenuations[i];
	}
	cout << "The total thickness of " << totalThickness << " could not be found in the attenuation list" << endl;
	return 0;
}


// Function to calculate the factorial of an integer recursively
int factorial(unsigned n)
{
	if (n == 0) return 1;
	return n*factorial(n-1);
}


// Function to transform an integer to a C++ string recursively
string intToString(int n)
{
	string number = "";
	if (n < 0) {	// Check if the integer is negative
		number += '-';
		number += intToString(-n);
	} else if (n < 10) {
		number += char(int('0')+n); // int('0') == 48 in ASCII
	} else {
		int mod = n % 10;
		int rest = (n-mod)/10;
		number += intToString(rest);
		number += char(int('0')+mod);
	}
	return number;
}
