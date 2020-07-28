#include "L1Trigger/L1TGlobal/interface/FractionalPrescalesVetosHelper.h"

using namespace l1t;

const FractionalPrescalesVetosHelper* FractionalPrescalesVetosHelper::readFromEventSetup(const L1TGlobalFractionalPrescalesVetos* es) {
  return new FractionalPrescalesVetosHelper(es);
}

FractionalPrescalesVetosHelper* FractionalPrescalesVetosHelper::readAndWriteFromEventSetup(const L1TGlobalFractionalPrescalesVetos* es) {
  FractionalPrescalesVetosHelper* x = new FractionalPrescalesVetosHelper(es);
  x->useCopy();
  return x;
}

FractionalPrescalesVetosHelper::FractionalPrescalesVetosHelper(L1TGlobalFractionalPrescalesVetos* w) {
  write_ = w;
  check_write();
  we_own_write_ = false;
  write_->version_ = VERSION_;
  read_ = write_;
}

FractionalPrescalesVetosHelper::FractionalPrescalesVetosHelper(const L1TGlobalFractionalPrescalesVetos* es) {
  read_ = es;
  write_ = nullptr;
}

void FractionalPrescalesVetosHelper::useCopy() {
  write_ = new L1TGlobalFractionalPrescalesVetos(*read_);
  we_own_write_ = true;
  read_ = write_;
}

FractionalPrescalesVetosHelper::~FractionalPrescalesVetosHelper() {
  if (we_own_write_ && write_)
    delete write_;
}
