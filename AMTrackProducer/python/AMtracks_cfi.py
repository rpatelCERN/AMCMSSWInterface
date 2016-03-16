import FWCore.ParameterSet.Config as cms

AMTrackProducer = cms.EDProducer('AMTrackProducer',
    inputTagStub=cms.InputTag('TTStubsFromPixelDigis', 'StubAccepted'),
    RoadsInputTag=cms.InputTag("AMRoadProducer", "Level1Roads"),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_ntuple*_*_*",
    )
)
