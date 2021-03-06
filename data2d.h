/*
 *  data2d.h
 *  cheetah
 *
 *  Created by Anton Barty on 4/9/10.
 *  Copyright 2010 all rights reserved.
 *
 *	You can modify this software under the terms of the GNU General Public License 
 *	as published by the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with this program.  If not, see <http://www.gnu.org/licenses/>. 
 */

#ifndef _data2d_h
#define _data2d_h

#define	tData2d	float

class cData2d {
	
public:
	cData2d();
	cData2d(long);
	~cData2d();
	
	void create(long);
	void create(long, long);
	void readHDF5(char *);
	void readHDF5(char *, char*);
	void writeHDF5(char *);

	
public:
	long		nx,ny,nn;
	tData2d		*data;
	
	
private:
	
	
};

#endif
