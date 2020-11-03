#include <iomanip>
#include <fstream>
#include <iostream>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "CondFormats/DataRecord/interface/L1TGlobalPrescalesVetosFractRcd.h"
#include "CondFormats/L1TObjects/interface/L1TGlobalPrescalesVetosFract.h"
#include "L1Trigger/L1TGlobal/interface/PrescalesVetosFractHelper.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"

class L1TPrescaleDumper : public edm::EDAnalyzer {
public:
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  explicit L1TPrescaleDumper(const edm::ParameterSet&) : edm::EDAnalyzer() {}
  ~L1TPrescaleDumper(void) override {}
};

void L1TPrescaleDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& evSetup) {

  edm::ESHandle<L1TGlobalPrescalesVetosFract> handle1;
  evSetup.get<L1TGlobalPrescalesVetosFractRcd>().get(handle1);

  const L1TGlobalPrescalesVetosFract* es = handle1.product();
  const l1t::PrescalesVetosFractHelper* m_prescale_helper = l1t::PrescalesVetosFractHelper::readFromEventSetup(es);

  std::vector<std::vector<double>> prescaleset = m_prescale_helper->prescaleTable();

  for(unsigned int i=0; i<prescaleset.size(); i++){

    const std::vector<double>& prescaleAlgoTrig = prescaleset.at(i);

    std::cout << "Index: " << i << std::endl;

    for(unsigned int j=0; j<prescaleAlgoTrig.size(); j++){
      std::cout << "L1T Algo #"<< j << ": "<< prescaleAlgoTrig.at(j) << std::endl;
    }
  }
}


#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"

DEFINE_FWK_MODULE(L1TPrescaleDumper);
