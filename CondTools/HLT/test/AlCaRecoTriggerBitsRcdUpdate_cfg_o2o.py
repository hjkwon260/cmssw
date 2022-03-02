# Config file template to write new/update AlCaRecoTriggerBits stored
# in AlCaRecoTriggerBitsRcd that is used to get selected HLT paths for
# the HLTHighLevel filter for AlCaReco production.
#
# Please understand that there are two IOVs involved:
# 1) One for the output tag. Here the usually used default is 1->inf,
#    changed by process.AlCaRecoTriggerBitsRcdUpdate.firstRunIOV
#    and process.AlCaRecoTriggerBitsRcdUpdate.lastRunIOV.
# 2) The IOV of the tag of the input AlCaRecoTriggerBitsRcd.
#    That is chosen by process.source.firstRun (but irrelevant if 
#    process.AlCaRecoTriggerBitsRcdUpdate.startEmpty = True)
#
# See also further comments below, especially the WARNING.
#
#  Author    : Marco Musich
#  Date      : Feb 2016
#  Modified  : Hyejin Kwon
#  Date      : Nov 2021

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing 

process = cms.Process("UPDATEDB")

options = VarParsing.VarParsing()
options.register( "inputDB", 
                  "frontier://FrontierProd/CMS_CONDITIONS",  #default value
                  VarParsing.VarParsing.multiplicity.singleton, 
                  VarParsing.VarParsing.varType.string,
                  "the input DB"
                  )

options.register( "inputTag", 
                  "AlCaRecoHLTpaths8e29_5e33_v7_prompt",  #default value
                  VarParsing.VarParsing.multiplicity.singleton, 
                  VarParsing.VarParsing.varType.string,
                  "the input tag"
                  )

options.register( "outputDB", 
                  "sqlite_file:AlCaRecoTriggerBits.db",  #default value
                  VarParsing.VarParsing.multiplicity.singleton, 
                  VarParsing.VarParsing.varType.string,
                  "the output DB"
                  )

options.register( "outputTag", 
                  "AlCaRecoTriggerBitsTag",  #default value
                  VarParsing.VarParsing.multiplicity.singleton, 
                  VarParsing.VarParsing.varType.string,
                  "the output tag"
                  )

options.register( "firstRun", 
                  1,  #default value
                  VarParsing.VarParsing.multiplicity.singleton, 
                  VarParsing.VarParsing.varType.int,
                  "the first run"
                  )

options.register( "hltKey", 
                  "/cdaq/special/commissioning2021/CRAFT/Cosmics/V4",  #default value
                  VarParsing.VarParsing.multiplicity.singleton, 
                  VarParsing.VarParsing.varType.string,
                  "the hlt key"
                  )

options.register( "keyToModify", 
                  "SiStripCalMinBias",  #default value
                  VarParsing.VarParsing.multiplicity.singleton, 
                  VarParsing.VarParsing.varType.string,
                  "the key to modify"
                  )

options.register('pathsToModify',
                 'HLT_*ZeroBias_part*_v*', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Comma-separated list of paths to be modified")
options.parseArguments()

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr = cms.untracked.PSet(enable = cms.untracked.bool(True))
process.MessageLogger.cout = cms.untracked.PSet(INFO = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(1)
    ))
process.MessageLogger.cout.enable = cms.untracked.bool(True)
process.MessageLogger.cout.threshold = cms.untracked.string('DEBUG')
process.MessageLogger.debugModules = cms.untracked.vstring('*')
# the module writing to DB
from CondTools.HLT.alCaRecoTriggerBitsRcdUpdate_cfi import alCaRecoTriggerBitsRcdUpdate as _alCaRecoTriggerBitsRcdUpdate
process.AlCaRecoTriggerBitsRcdUpdate = _alCaRecoTriggerBitsRcdUpdate.clone()
# The IOV that you want to write out, defaut is 1 to -1/inf. 
process.AlCaRecoTriggerBitsRcdUpdate.firstRunIOV = options.firstRun # docu see...
#process.AlCaRecoTriggerBitsRcdUpdate.lastRunIOV = -1 # ...cfi
# If you want to start from scratch, comment the next line:
process.AlCaRecoTriggerBitsRcdUpdate.startEmpty = False

# add paths if hlt key has 'special', remove those for the other cases
runMode = options.hltKey.split('/')[2]
if(runMode=='special'):
  print(runMode,'mode: adding',options.pathsToModify,'to',options.keyToModify,'if not exist')
  process.AlCaRecoTriggerBitsRcdUpdate.pathsToAdd = [
        cms.PSet(listName = cms.string(options.keyToModify),
                 hltPaths = cms.vstring(options.pathsToModify.split(','))
                 )
  ]
else: 
  print('not special mode, removing',options.pathsToModify,'from',options.keyToModify,'if exist') 
  process.AlCaRecoTriggerBitsRcdUpdate.pathsToRemove = [
        cms.PSet(listName = cms.string(options.keyToModify),
                 hltPaths = cms.vstring(options.pathsToModify.split(','))
                 )
  ]

process.source = cms.Source("EmptySource",
                            firstRun = cms.untracked.uint32(options.firstRun) # runnumber
                            )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("CondCore.CondDB.CondDB_cfi")

# DB input service: 
process.CondDB.connect = options.inputDB
process.dbInput = cms.ESSource("PoolDBESSource",
                               process.CondDB,
                               toGet = cms.VPSet(cms.PSet(record = cms.string('AlCaRecoTriggerBitsRcd'),
                                                          tag = cms.string(options.inputTag)
                                                          )
                                                 )
                               )

# DB output service:
process.CondDB.connect = options.outputDB
process.PoolDBOutputService = cms.Service("PoolDBOutputService",
                                          process.CondDB,
                                          timetype = cms.untracked.string('runnumber'),
                                          toPut = cms.VPSet(cms.PSet(record = cms.string('AlCaRecoTriggerBitsRcd'),
                                                                     tag = cms.string(options.outputTag) # choose output tag you want
                                                                     )
                                                            )
                                          )

# Put module in path:
process.p = cms.Path(process.AlCaRecoTriggerBitsRcdUpdate)


