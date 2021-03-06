#ifndef Pds_EvrConfigType_hh
#define Pds_EvrConfigType_hh

#include "pdsdata/xtc/TypeId.hh"
#include "pdsdata/evr/ConfigV7.hh"
#include "pdsdata/evr/DataV3.hh"

typedef Pds::EvrData::ConfigV7  EvrConfigType;
typedef Pds::EvrData::DataV3    EvrDataType;

static Pds::TypeId _evrConfigType(Pds::TypeId::Id_EvrConfig,
          EvrConfigType::Version);

#endif
