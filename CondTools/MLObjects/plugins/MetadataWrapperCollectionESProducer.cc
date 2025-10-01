#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetupRecordIntervalFinder.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/ESInputTag.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"

#include "CondFormats/MLObjects/interface/Metadata.h"
#include "CondFormats/DataRecord/interface/MetadataRcd.h"
#include "CondFormats/MLObjects/interface/MetadataWrapper.h"
#include "CondFormats/DataRecord/interface/MetadataWrapperRcd.h"

#include "nlohmann/json.hpp"
#include <fstream>
#include <iostream>

class MetadataWrapperCollectionESProducer : public edm::ESProducer {
public:
  MetadataWrapperCollectionESProducer(const edm::ParameterSet&);
  std::unique_ptr<MetadataWrapperCollection> produce(const MetadataWrapperRcd&);

private:
  edm::ESGetToken<MetadataCollection, MetadataRcd> token_;
  nlohmann::json json_;
};

MetadataWrapperCollectionESProducer::MetadataWrapperCollectionESProducer(const edm::ParameterSet& iConfig) {

  std::string label = iConfig.getParameter<std::string>("label");
  std::string jsonFile = iConfig.getParameter<std::string>("jsonFile");

  auto cc = setWhatProduced(this, label);

  token_ = cc.consumesFrom<MetadataCollection, MetadataRcd>(edm::ESInputTag{"", label});

  // Load JSON file once at construction
  std::ifstream in(jsonFile);
  if (!in.is_open()) {
    throw cms::Exception("FileNotFound") << "Could not open JSON file: " << jsonFile;
  }
  in >> json_;

}

std::unique_ptr<MetadataWrapperCollection> MetadataWrapperCollectionESProducer::produce(const MetadataWrapperRcd& iRecord) {

  MetadataCollection const&  cond = iRecord.get(token_);

  auto wrapper = std::make_unique<MetadataWrapperCollection>();

  for (auto const& m : cond.models()) {

      std::string model_path = "";
      std::string preproc_path = "";

      for (auto const& entry : json_) {
          if (entry["model_name"] == m.model_name() &&
          entry["version"] == std::to_string(m.version())) {
          model_path = entry["model_path"].get<std::string>();
          preproc_path = entry["preprocessing_path"].get<std::string>();
          break;
        }
      }

      wrapper->add(MetadataWrapper(m.model_name(), std::to_string(m.version()), model_path, preproc_path));
  }


  return wrapper;

}

DEFINE_FWK_EVENTSETUP_MODULE(MetadataWrapperCollectionESProducer);