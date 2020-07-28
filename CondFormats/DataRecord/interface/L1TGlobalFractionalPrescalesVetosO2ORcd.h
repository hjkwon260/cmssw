// L1TGlobalPrescalesVetosRcd
// Description: Record for L1TGlobalPrescalesVetos
//
// automatically generate by make_records.pl
//
#ifndef CondFormatsDataRecord_L1TGlobalFractionalPrescalesVetosO2O_h
#define CondFormatsDataRecord_L1TGlobalFractionalPrescalesVetosO2O_h

#include "FWCore/Framework/interface/DependentRecordImplementation.h"
#include "CondFormats/DataRecord/interface/L1TriggerKeyListExtRcd.h"
#include "CondFormats/DataRecord/interface/L1TriggerKeyExtRcd.h"
#include "CondFormats/DataRecord/interface/L1TGlobalFractionalPrescalesVetosRcd.h"
class L1TGlobalFractionalPrescalesVetosO2ORcd
    : public edm::eventsetup::DependentRecordImplementation<
          L1TGlobalFractionalPrescalesVetosO2ORcd,
          boost::mpl::vector<L1TriggerKeyListExtRcd, L1TriggerKeyExtRcd, L1TGlobalFractionalPrescalesVetosRcd> > {};

#endif
