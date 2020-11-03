import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("tester")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cout.placeholder = cms.untracked.bool(False)
process.MessageLogger.cout.threshold = cms.untracked.string('DEBUG')
process.MessageLogger.debugModules = cms.untracked.vstring('*')

process.source = cms.Source("EmptySource", firstRun = cms.untracked.uint32(337857))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

from CondCore.CondDB.CondDB_cfi import CondDB
CondDB.connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS')

process.l1conddb = cms.ESSource("PoolDBESSource",
       CondDB,
       toGet   = cms.VPSet(
            cms.PSet(
                 record = cms.string('L1TGlobalPrescalesVetosFractRcd'),
                 tag = cms.string("L1TGlobalPrescalesVetosFract_Stage2v1_hlt")
            )
       )
)

process.l1pu = cms.EDAnalyzer("L1TPrescaleDumper")

process.p = cms.Path(process.l1pu)

