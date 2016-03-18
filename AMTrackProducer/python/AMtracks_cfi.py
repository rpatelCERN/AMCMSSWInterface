import FWCore.ParameterSet.Config as cms

AMTrackProducer = cms.EDProducer('AMTrackProducer',
    inputTagStub=cms.InputTag('TTStubsFromPixelDigis', 'StubAccepted'),
    # RoadsInputTag=cms.InputTag("AMRoadProducer", "Level1Roads"),
    RoadsInputTag= cms.InputTag("MergeTCOutput", "AML1TCs"),
    ConstantsDir=cms.FileInPath("AMCMSSWInterface/AMTrackProducer/python/ConstantsProduction/PreEstimate_Transverse/matrixVD_2016.txt"),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_ntuple*_*_*",
    )
)
