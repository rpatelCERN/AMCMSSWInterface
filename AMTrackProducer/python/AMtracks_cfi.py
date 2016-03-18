import FWCore.ParameterSet.Config as cms

AMTrackProducer = cms.EDProducer('AMTrackProducer',
    inputTagStub=cms.InputTag('TTStubsFromPixelDigis', 'StubAccepted'),
    RoadsInputTag= cms.InputTag("MergeTCOutput", "AML1TCs"),
    #RoadsInputTag=cms.InputTag("AMRoadProducer", "Level1Roads"),
    ConstantsDir=cms.string("/fdata/hepx/store/user/rish/AMSIMULATION/Forked/CMSSW_6_2_0_SLHC25_patch3/src/AMCMSSWInterface/AMTrackProducer/python/ConstantsProduction/"),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_ntuple*_*_*",
    )
)
