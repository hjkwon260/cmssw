#ifndef FractionalPRESCALESVETOSHELPERS_H__
#define FractionalPRESCALESVETOSHELPERS_H__

#include <cassert>
#include "CondFormats/L1TObjects/interface/L1TGlobalFractionalPrescalesVetos.h"

// If you want to create a new object that you can read and write, use this constructor:
//
//   l1t::PrescalesVetosHelper x(new L1TPrescalesVetors());
//
// If you wish to read the table from the EventSetup, and will only read, use this:
//
//   const PrescalesVetosHelper * x = PrescalesVetosHelper::readFromEventSetup(...)
//   //...
//   delete x;
//
// If you wish to read the table from the EventSetup, but then be able to edit the values locally, use this:
//
//   PrescalesVetorsHelper * x = PrescalesVetosHelper::readAndWriteFromEventSetup(...)
//   //...
///  delete x;
//
// but there's a performance penalty as a copy is made.

//
// This class does not take over responsibility for deleting the pointers it is
// initialized with.  That is responsibility of the calling code.
//

namespace l1t {

  class FractionalPrescalesVetosHelper {
  public:
    enum { VERSION_ = 1 };

    ~FractionalPrescalesVetosHelper();

    //ctor if creating a new table (e.g. from XML or python file)
    FractionalPrescalesVetosHelper(L1TGlobalFractionalPrescalesVetos* w);
    //create for reading only, from the EventSetup:
    static const FractionalPrescalesVetosHelper* readFromEventSetup(const L1TGlobalFractionalPrescalesVetos* es);
    // create for reading and writing, starting from the EventSetup:
    static FractionalPrescalesVetosHelper* readAndWriteFromEventSetup(const L1TGlobalFractionalPrescalesVetos* es);

    int bxMaskDefault() const { return read_->bxmask_default_; };
    void setBxMaskDefault(int value) {
      check_write();
      write_->bxmask_default_ = value;
    };

    inline const std::vector<std::vector<double> >& prescaleTable() const { return read_->prescale_table_; };
    void setPrescaleFactorTable(std::vector<std::vector<double> > value) {
      check_write();
      write_->prescale_table_ = value;
    };
    inline const std::vector<int>& triggerMaskVeto() const { return read_->veto_; };
    void setTriggerMaskVeto(std::vector<int> value) {
      check_write();
      write_->veto_ = value;
    };

    inline const std::map<int, std::vector<int> >& triggerAlgoBxMask() const { return read_->bxmask_map_; };
    void setTriggerAlgoBxMask(std::map<int, std::vector<int> > value) {
      check_write();
      write_->bxmask_map_ = value;
    };

    // access to underlying pointers, mainly for ESProducer:
    const L1TGlobalFractionalPrescalesVetos* getReadInstance() const { return read_; }
    L1TGlobalFractionalPrescalesVetos* getWriteInstance() { return write_; }

  private:
    FractionalPrescalesVetosHelper(const L1TGlobalFractionalPrescalesVetos* es);
    void useCopy();
    void check_write() { assert(write_); }
    // separating read from write allows for a high-performance read-only mode (as no copy is made):
    const L1TGlobalFractionalPrescalesVetos* read_;  // when reading/getting, use this.
    L1TGlobalFractionalPrescalesVetos* write_;       // when writing/setting, use this.
    bool we_own_write_;
  };

}  // namespace l1t
#endif
