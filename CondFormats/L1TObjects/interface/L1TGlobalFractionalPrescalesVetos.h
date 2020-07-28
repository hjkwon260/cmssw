//  L1TGlobalPrescalesVetos
//
//  Table containing the entire set of prescales and masks for each L1T algorithm bit
//

#ifndef L1TGlobalFractionalPrescalesVetos_h
#define L1TGlobalFractionalPrescalesVetos_h

#include <vector>

#include "CondFormats/Serialization/interface/Serializable.h"

class L1TGlobalFractionalPrescalesVetos {
public:
  L1TGlobalFractionalPrescalesVetos() {
    version_ = 0;
    bxmask_default_ = 0;
  }

  unsigned int version_;
  std::vector<std::vector<double> > prescale_table_;
  int bxmask_default_;
  std::map<int, std::vector<int> > bxmask_map_;
  std::vector<int> veto_;
  std::vector<int> exp_ints_;
  std::vector<double> exp_doubles_;

  COND_SERIALIZABLE;
};

#endif
