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
import fnmatch, subprocess, time, re

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
                  "/cdaq/special/PilotBeamTest2021/Collisions/V55",  #default value
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
                 'HLT_ZeroBias_part*_v*', #default value
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

start = time.time()

print('HLT menu:', options.hltKey)

cmd = "hltGetConfiguration adg:%s" %options.hltKey
p = subprocess.Popen( cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
out, err = p.communicate()

# if hltGetConfiguration fails
if p.returncode != 0:
    raise Exception(out.decode('utf-8'))

# get part where pds:paths are defined
is_block = False
lines=''
for line in out.decode('utf-8').split('\n'):
  if is_block and 'process.' in line:
    break  
  if 'process.datasets' in line:
    is_block = True
  if is_block:
    lines += line+"\n"

exec(lines)

pds = process.datasets._Parameterizable__parameterNames
for pd in pds:
  paths = sorted( path for path in process.datasets.__dict__[pd] )
  if (re.search("^ZeroBias\d{1,2}", pd)): # only search for 'ZeroBiasX'
    pattern = 'HLT_*ZeroBias*part*'
    matching = fnmatch.filter(paths, pattern)
    if matching:
      print(pd, ':', *matching)

# add partitioned paths in ZeroBiasX if exist, remove those for the other cases
if 'matching' in locals():
  print('adding',options.pathsToModify,'to',options.keyToModify,'if not exist')
  process.AlCaRecoTriggerBitsRcdUpdate.pathsToAdd = [
        cms.PSet(listName = cms.string(options.keyToModify),
                 hltPaths = cms.vstring(options.pathsToModify.split(','))
                 )
  ]
else: 
  print('removing',options.pathsToModify,'from',options.keyToModify,'if exist') 
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

end = time.time()
print('Elapsed time for O2O: ', end - start, 'sec')
