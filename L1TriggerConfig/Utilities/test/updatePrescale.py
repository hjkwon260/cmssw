import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')

options.register('iov',
                  1, # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.int,         # string, int, or float
                  "iov")
options.parseArguments()

process = cms.Process("tester")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cout.placeholder = cms.untracked.bool(False)
process.MessageLogger.cout.threshold = cms.untracked.string('DEBUG')
process.MessageLogger.debugModules = cms.untracked.vstring('*')

process.source = cms.Source("EmptySource", firstRun = cms.untracked.uint32(options.iov))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

from CondCore.CondDB.CondDB_cfi import CondDB
CondDB.connect = cms.string('sqlite:l1config.db')

process.l1conddb = cms.ESSource("PoolDBESSource",
       CondDB,
       toGet   = cms.VPSet(
            cms.PSet(
                 record = cms.string('L1TGlobalPrescalesVetosRcd'),
                 tag = cms.string("L1TGlobalPrescalesVetos_Stage2v0_hlt")
            )
       )
)

process.l1pu = cms.EDAnalyzer("L1TPrescaleUpdater")

from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
outputDB = cms.Service("PoolDBOutputService",
                       CondDBSetup,
                       connect = cms.string('sqlite_file:l1configTweak.db'),
                       toPut   = cms.VPSet(
                           cms.PSet(
                               record = cms.string('L1TGlobalPrescalesVetosFractRcd'),
                               tag = cms.string("L1TGlobalPrescalesVetosFract_Stage2v0_hlt")
                           )
                       )
)
outputDB.DBParameters.authenticationPath = '.'
process.add_(outputDB)

process.p = cms.Path(process.l1pu)

