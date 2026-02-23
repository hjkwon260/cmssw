import FWCore.ParameterSet.Config as cms

from RecoBTag.FeatureTools.pfDeepBoostedJetTagInfos_cfi import pfDeepBoostedJetTagInfos
from RecoBTag.ONNXRuntime.boostedJetONNXJetTagsProducer_cfi import boostedJetONNXJetTagsProducer
process = cms.Process("TESTML")

# Minimal services
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Empty input source (no events needed for metadata test)
process.source = cms.Source("EmptySource")
process.source.firstRun = cms.untracked.uint32(123456)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

# CondDB connection
process.load("CondCore.CondDB.CondDB_cfi")
process.CondDB.connect = 'sqlite_file:MLMetadata.db'
process.PoolDBESSource = cms.ESSource("PoolDBESSource",
    process.CondDB,
    toGet = cms.VPSet(
        cms.PSet(
            record = cms.string("MetadataRcd"),
            tag    = cms.string("MLMetadata_BTV_UParT"),
            label  = cms.untracked.string("BTV_UParT") 
        )      
    )
)

process.mlModelPathResolver = cms.ESProducer(
    "MetadataWrapperESProducer",
    label = cms.string("BTV_UParT"),
    jsonFile = cms.string("RecoBTag/Combined/RobustParTAK4/PUPPI/V00/modelfile/UParTAK4.json")  # local JSON file        
)


process.pfParticleNetJetTags = boostedJetONNXJetTagsProducer.clone(
    src = 'pfParticleNetTagInfos',
    preprocess_json = 'RecoBTag/Combined/data/ParticleNetAK8/General/V01/preprocess.json',
    model_path = 'RecoBTag/Combined/data/ParticleNetAK8/General/V01/modelfile/model.onnx',
    flav_names = ["probTbcq",  "probTbqq",  "probTbc",   "probTbq",  "probTbel", "probTbmu", "probTbta",
                  "probWcq",   "probWqq",   "probZbb",   "probZcc",  "probZqq",  "probHbb", "probHcc",
                  "probHqqqq", "probQCDbb", "probQCDcc", "probQCDb", "probQCDc", "probQCDothers"],
    label  = cms.untracked.string("BTV_UParT")
)

# --- Path ---
process.p = cms.Path(process.pfParticleNetJetTags)

