#include <iomanip>
#include <fstream>
#include <iostream>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "CondFormats/DataRecord/interface/L1TGlobalPrescalesVetosRcd.h"
#include "CondFormats/DataRecord/interface/L1TGlobalPrescalesVetosFractRcd.h"
#include "CondFormats/L1TObjects/interface/L1TGlobalPrescalesVetos.h"
#include "CondFormats/L1TObjects/interface/L1TGlobalPrescalesVetosFract.h"
#include "L1Trigger/L1TGlobal/interface/PrescalesVetosHelper.h"
#include "L1Trigger/L1TGlobal/interface/PrescalesVetosFractHelper.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"

class L1TPrescaleUpdater : public edm::EDAnalyzer {
public:
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  explicit L1TPrescaleUpdater(const edm::ParameterSet&) : edm::EDAnalyzer() {}
  ~L1TPrescaleUpdater(void) override {}
};

void L1TPrescaleUpdater::analyze(const edm::Event& iEvent, const edm::EventSetup& evSetup) {

  edm::ESHandle<L1TGlobalPrescalesVetos> handle1;
  evSetup.get<L1TGlobalPrescalesVetosRcd>().get(handle1);

  const L1TGlobalPrescalesVetos* es = handle1.product();
  const l1t::PrescalesVetosHelper* m_prescale_helper = l1t::PrescalesVetosHelper::readFromEventSetup(es);

  //convert int prescale to double
  std::vector<std::vector<int>> prescales = m_prescale_helper->prescaleTable();
  std::vector<std::vector<double>> prescales_fract;
  prescales_fract.reserve(prescales.size());

  for (auto&& v : prescales) prescales_fract.emplace_back(std::begin(v), std::end(v));  

  l1t::PrescalesVetosFractHelper data_(new L1TGlobalPrescalesVetosFract());

  data_.setBxMaskDefault(m_prescale_helper->bxMaskDefault());
  data_.setPrescaleFactorTable(prescales_fract);
  data_.setTriggerMaskVeto(m_prescale_helper->triggerMaskVeto());
  data_.setTriggerAlgoBxMask(m_prescale_helper->triggerAlgoBxMask());
  
  auto payload = std::make_unique<const L1TGlobalPrescalesVetosFract>(*data_.getWriteInstance());

  edm::Service<cond::service::PoolDBOutputService> poolDb;
  if (poolDb.isAvailable()) {
    cond::Time_t currentTime = poolDb->currentTime();
    poolDb->writeOne(payload.get(), currentTime, "L1TGlobalPrescalesVetosFractRcd");

  }
}


#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"

DEFINE_FWK_MODULE(L1TPrescaleUpdater);
