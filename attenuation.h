/*
 *  attenuation.h
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

#ifndef _attenuation_h
#define _attenuation_h

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

/*
 *	Function prototypes for calculating attenuation
 */

int getSiThickness(unsigned & totalThickness, unsigned nFilters, unsigned filterThicknesses[]);
double getAttenuation(unsigned totalThickness, unsigned nThicknesses, unsigned possibleThicknesses[], double possibleAttenuations[]);
int factorial(unsigned n);
std::string intToString(int n);
template <class T>
bool fromString(T & t, const std::string & s, std::ios_base & (*f)(std::ios_base &));

/*
 *	Implementations of templatized functions used for calculating attenuation
 */

// Templatized function to transform a number to a C++ string
template <class T>
bool fromString(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&))
{
	std::istringstream iss(s);
	return !(iss >> f >> t).fail();
}

#endif
