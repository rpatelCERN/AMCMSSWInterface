import FWCore.ParameterSet.Config as cms

AMRoadProducer = cms.EDProducer('AMRoadProducer',
    inputTagStub=cms.InputTag('TTStubsFromPixelDigis', 'StubAccepted'),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_ntuple*_*_*",
    )
)
