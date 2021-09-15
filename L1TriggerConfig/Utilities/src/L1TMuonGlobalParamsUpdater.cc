#include <iomanip>
#include <iostream>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "CondFormats/DataRecord/interface/L1TMuonGlobalParamsO2ORcd.h"
#include "CondFormats/DataRecord/interface/L1TMuonGlobalParamsRcd.h"
#include "CondFormats/L1TObjects/interface/L1TMuonGlobalParams.h"
#include "L1Trigger/L1TMuon/interface/L1TMuonGlobalParamsHelper.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"

class L1TMuonGlobalParamsUpdater : public edm::EDAnalyzer {
public:
    void analyze(const edm::Event&, const edm::EventSetup&) override;

    explicit L1TMuonGlobalParamsUpdater(const edm::ParameterSet &pset) : edm::EDAnalyzer(){
    }
    ~L1TMuonGlobalParamsUpdater(void) override{}
};

void L1TMuonGlobalParamsUpdater::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
    edm::ESHandle<L1TMuonGlobalParams> handle1;
    iSetup.get<L1TMuonGlobalParamsRcd>().get(handle1);

    const L1TMuonGlobalParams* es = handle1.product();

    std::unique_ptr<L1TMuonGlobalParamsHelper> m_param_helper = std::make_unique<L1TMuonGlobalParamsHelper>();

    m_param_helper = std::make_unique<L1TMuonGlobalParamsHelper>(*es);

    m_param_helper->setFwVersion(0x6010000);

    auto retval = std::make_unique<const L1TMuonGlobalParams>(cast_to_L1TMuonGlobalParams((L1TMuonGlobalParams_PUBLIC)*m_param_helper));

    edm::Service<cond::service::PoolDBOutputService> poolDb;
    if( poolDb.isAvailable() ){
        cond::Time_t firstSinceTime = poolDb->beginOfTime();
        poolDb->writeOne(retval.get(),firstSinceTime,"L1TMuonGlobalParamsRcd");
    }
}

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"

DEFINE_FWK_MODULE(L1TMuonGlobalParamsUpdater);

