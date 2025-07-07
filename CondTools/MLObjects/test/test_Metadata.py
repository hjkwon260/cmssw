import FWCore.ParameterSet.Config as cms

process = cms.Process("MDWriter")
process.load("CondCore.CondDB.CondDB_cfi")
process.CondDB.connect = 'sqlite_file:MLMetadata.db'

process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(1),
)

process.source = cms.Source("EmptySource")
process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    process.CondDB,
    toPut = cms.VPSet(
        cms.PSet(
            record = cms.string("MetadataRcd"),
            tag = cms.string("MLMetadata_Tag"),
            label = cms.string(""),
        ),
    )
)

process.dbCreator = cms.EDAnalyzer("MetadataWriter",
    value = cms.int32(1993),
    info = cms.string("test ML Metadata DB")
)

process.p = cms.Path(
    process.dbCreator
)
