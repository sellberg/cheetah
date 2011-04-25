#ifndef _commonmode_h
#define _commonmode_h
/*
 *  commonmode.h
 *  cheetah
 *
 *  Created by Anton Barty on 6/2/11.
 *  Copyright 2011 CFEL. All rights reserved.
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

#include "setup.h"
#include "worker.h"

/*
 *	Function prototypes
 */
void cmModuleSubtract(tThreadInfo*, cGlobal*);
void cmSubModuleSubtract(tThreadInfo*, cGlobal*);


#endif
