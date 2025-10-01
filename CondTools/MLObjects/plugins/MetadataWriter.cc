#include <iostream>
#include <fstream>
#include <vector>
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

  const std::vector<edm::ParameterSet> models_;  
};

MetadataWriter::MetadataWriter(const edm::ParameterSet& iConfig)
    : models_(iConfig.getParameter<std::vector<edm::ParameterSet>>("models")) {}

void MetadataWriter::beginJob() {

    MetadataCollection coll;

    for (auto const& ps : models_) {
        int version = ps.getParameter<int>("version");
        std::string model_name = ps.getParameter<std::string>("model_name");
        coll.add_model(Metadata(model_name, version));
    }

  edm::Service<cond::service::PoolDBOutputService> pool;
  if (pool.isAvailable())
        pool->writeOneIOV(coll, pool->currentTime(), "MetadataRcd");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MetadataWriter);  