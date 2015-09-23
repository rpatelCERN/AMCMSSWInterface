import FWCore.ParameterSet.Config as cms
from SLHCUpgradeSimulations.L1TrackTrigger.l1TkMuonsExt_cfi import l1TkMuonsExt

l1TkMuonsExtCSC = l1TkMuonsExt.clone(
                 L1MuonsInputTag = cms.InputTag("l1extraMuExtended", "csc"),
                 L1TrackInputTag = cms.InputTag("AMTrackProducer", "Level1TTTracks")
		)
