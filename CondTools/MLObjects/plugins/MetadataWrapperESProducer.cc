/**
 * Author: H. Kwon
 */

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
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"

class MetadataWrapperESProducer : public edm::ESProducer {
public:
  MetadataWrapperESProducer(const edm::ParameterSet&);
  std::unique_ptr<MetadataWrapper> produce(const MetadataWrapperRcd&);

private:
  edm::ESGetToken<Metadata, MetadataRcd> token_;
  nlohmann::json json_;
};

MetadataWrapperESProducer::MetadataWrapperESProducer(const edm::ParameterSet& iConfig)
{
  std::string label = iConfig.getParameter<std::string>("label");
  std::string jsonFile = iConfig.getParameter<std::string>("jsonFile");

  auto cc = setWhatProduced(this, label);
  token_ = cc.consumesFrom<Metadata, MetadataRcd>(edm::ESInputTag{"", label});

  // Load JSON file once at construction
  std::ifstream in(jsonFile);
  if (!in.is_open()) {
    throw cms::Exception("FileNotFound") << "Could not open JSON file: " << jsonFile;
  }
  in >> json_;
}

std::unique_ptr<MetadataWrapper> MetadataWrapperESProducer::produce(const MetadataWrapperRcd& iRecord) {
  const Metadata& cond = iRecord.get(token_);
  auto wrapper = std::make_unique<MetadataWrapper>();

  std::string model_path;
  std::string preproc_path;
  std::vector<std::string> input_features;
  std::vector<std::string> flav_names;

  for (const auto& entry : json_) {
      if (entry["model_name"] == cond.model_name() &&
          entry["version"] == std::to_string(cond.version()) &&
          entry["hash"] == cond.hash()) {
          model_path = entry["model_path"].get<std::string>();
          preproc_path = entry["preprocessing_path"].get<std::string>();
          if (entry.contains("input_features")) {
              input_features = entry["input_features"].get<std::vector<std::string>>();
          }
          if (entry.contains("flav_names")) {
              flav_names = entry["flav_names"].get<std::vector<std::string>>();
          }
          break;
      }
  }

  MetadataWrapper mw(cond.model_name(), std::to_string(cond.version()), cond.hash(), model_path, preproc_path);
  mw.set_input_features(input_features);
  mw.set_flav_names(flav_names);

  if (!model_path.empty()) {
      try {
          auto runtime = std::make_shared<cms::Ort::ONNXRuntime>(model_path);
          mw.setOnnxRuntime(runtime);
      } catch (const std::exception& e) {
          std::cout << "MetadataWrapperESProducer: Failed to load ONNX model " << cond.model_name()
                                                                << " v" << cond.version() << " h" << cond.hash() << ": " << e.what() << std::endl;
      }
  }

  return std::make_unique<MetadataWrapper>(mw);
}

DEFINE_FWK_EVENTSETUP_MODULE(MetadataWrapperESProducer);
