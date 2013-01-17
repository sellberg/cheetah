/*
 * CsPadConfigType.hh
 *
 *  Created on: Jun 28, 2010
 */

#ifndef CSPADCONFIGTYPE_HH_
#define CSPADCONFIGTYPE_HH_

#include "pdsdata/xtc/TypeId.hh"
#include "pdsdata/cspad/ConfigV4.hh"

typedef Pds::CsPad::ConfigV4 CsPadConfigType;

static Pds::TypeId _CsPadConfigType(Pds::TypeId::Id_CspadConfig,
                                    CsPadConfigType::Version);


#endif /* CSPADCONFIGTYPE_HH_ */
