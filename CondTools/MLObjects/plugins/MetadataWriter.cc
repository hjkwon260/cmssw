// system include files
#include <iostream>
#include <fstream>
#include <vector>

// user include files
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "CondFormats/DataRecord/interface/MetadataRcd.h"
#include "CondFormats/MLObjects/interface/Metadata.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"

class MetadataWriter : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:

  explicit MetadataWriter(const edm::ParameterSet&);
  ~MetadataWriter() override {}

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override {}

  // const int value_;
  // const std::string info_;
  const std::vector<edm::ParameterSet> models_;  // VPSet for multiple models
};

MetadataWriter::MetadataWriter(const edm::ParameterSet& iConfig)
    // : value_(iConfig.getParameter<int>("value")),
    //   info_(iConfig.getParameter<std::string>("info")){}
    : models_(iConfig.getParameter<std::vector<edm::ParameterSet>>("models")) {}

void MetadataWriter::beginJob() {

  // Metadata MD(value_, info_);
    MetadataCollection coll;
    // Fill collection from Python VPSet
    for (auto const& ps : models_) {
        int version = ps.getParameter<int>("version");
        std::string name = ps.getParameter<std::string>("name");
        coll.addModel(Metadata(version, name));
    }

  edm::Service<cond::service::PoolDBOutputService> pool;
  if (pool.isAvailable())
    // pool->writeOneIOV(MD, pool->currentTime(), "MetadataRcd");
        pool->writeOneIOV(coll, pool->currentTime(), "MetadataRcd");
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MetadataWriter);  