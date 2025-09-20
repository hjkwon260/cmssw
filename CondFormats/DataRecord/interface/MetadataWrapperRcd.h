// #ifndef DataRecord_MetadataWrapperRcd_h
// #define DataRecord_MetadataWrapperRcd_h

// #include "FWCore/Framework/interface/EventSetupRecordImplementation.h"

// class MetadataWrapperRcd : public edm::eventsetup::EventSetupRecordImplementation<MetadataWrapperRcd> {};

// #endif

#ifndef DataRecord_MetadataWrapperRcd_h
#define DataRecord_MetadataWrapperRcd_h

#include "FWCore/Framework/interface/DependentRecordImplementation.h"
#include "CondFormats/DataRecord/interface/MetadataRcd.h"
#include "FWCore/Utilities/interface/mplVector.h"

// Declare this record as dependent on MetadataRcd
class MetadataWrapperRcd : public edm::eventsetup::DependentRecordImplementation<
    MetadataWrapperRcd, edm::mpl::Vector<MetadataRcd> > {};

#endif