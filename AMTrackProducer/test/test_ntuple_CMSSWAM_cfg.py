import FWCore.ParameterSet.Config as cms
import sys
process = cms.Process("TEST2")
runOnMC = True
#use TreeMaker Options
from AMCMSSWInterface.AMTrackProducer.CommandLineParams import CommandLineParams

## MessageLogger
parameters = CommandLineParams()
inputFiles = parameters.value("inputFiles",
'file:/fdata/hepx/store/user/rish/AMSIMULATION/CMSSW_6_2_0_SLHC27/src/AMCMSSWInterface/AMTrackProducer/test/AMTC_output.root',
#'file:/fdata/hepx/store/user/rish/AMSIMULATION/Forked2/CMSSW_6_2_0_SLHC25_patch3/src/AMCMSSWInterface/AMTrackProducer/test/example_L1Roads.root',
#'file:/fdata/hepx/store/user/rish/CMSSW_6_2_0_SLHC25_patch3/src/Muons/SingleMuonNoPU_tt27_16.root'
#'file:/fdata/hepx/store/user/rish/Muons/TTIMuons11.root'
#'file:/fdata/hepx/store/user/rish/CMSSW_6_2_0_SLHC25_patch3/src/SingleMuonPU_tt27_%s.root'
#'file:/fdata/hepx/store/user/rish/CMSSW_6_2_0_SLHC25_patch3/src/SingleMuonPU_tt27_%s.root'
)
outputFile=parameters.value("outputFile","test_ntuple.root")
mode=parameters.value("mode","AM")
maxEvents=parameters.value("maxEvents", -1)
isoCuts=parameters.value("isoCuts",0.1)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1


## Source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(inputFiles)
)
## Maximal Number of Events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(maxEvents) )

## Geometry and Global Tags
process.load('Configuration.Geometry.GeometryExtended2023TTIReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023TTI_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')

## Write the TTree
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(outputFile)
)
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PH2_1K_FB_V3::All', '')

process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.Geometry.GeometryExtended2023TTIReco_cff')
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')
process.load('SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff')
process.load("AMCMSSWInterface.AMTrackProducer.AMTrackingSequence_cff")
#process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TkMuonSequence_cfi")
'''
if mode=="AM":
	process.load("AMCMSSWInterface.AMTrackProducer.L1TkMuonSequence_cfi")
else:
	process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TkMuonSequence_cfi")
'''
process.p = cms.Path(process.AMTracks)

#import the producers
#process.TTTrackAssociatorForAM=process.TTTrackAssociatorFromPixelDigis.clone()
#if mode == "AM":
#	process.TTTrackAssociatorForAM.TTTracks=cms.VInputTag( cms.InputTag("AMTrackProducer", "Level1TTTracks"))
#else:
#	process.TTTrackAssociatorForAM.TTTracks=cms.VInputTag( cms.InputTag("TTTracksFromPixelDigis", "Level1TTTracks"))




'''
process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TkEmParticleProducer_cfi")
process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TkElectronTrackProducer_cfi")
if mode == "AM":
	process.L1TkPhotons.L1TrackInputTag=cms.InputTag("AMTrackProducer","Level1TTTracks")

process.pL1TkPhotons = cms.Path( process.L1TkPhotons )
process.L1TkPhotonsIsolComp = process.L1TkPhotons.clone()
process.L1TkPhotonsIsolComp.IsoCut = cms.double( 0.10)
process.pL1TkPhotonsIsolTest = cms.Path( process.L1TkPhotonsIsolComp )
if mode == "AM":
	process.L1TkElectrons.L1TrackInputTag=cms.InputTag("AMTrackProducer","Level1TTTracks")

process.pElectrons = cms.Path( process.L1TkElectrons )

process.L1TkIsoElectrons = process.L1TkElectrons.clone()
process.L1TkIsoElectrons.IsoCut = cms.double(0.1)
process.pL1TkIsoElectrons=cms.Path(process.L1TkIsoElectrons)

process.L1TkElectronsLoose = process.L1TkElectrons.clone()
process.L1TkElectronsLoose.TrackEGammaDeltaPhi = cms.vdouble(0.07, 0.0, 0.0)
process.L1TkElectronsLoose.TrackEGammaDeltaR = cms.vdouble(0.12, 0.0, 0.0)
process.L1TkElectronsLoose.TrackMinPt = cms.double( 3.0 )
process.pElectronsLoose = cms.Path( process.L1TkElectronsLoose)
AMTrackInputTag = cms.InputTag("AMTrackProducer", "Level1TTTracks")
'''
#if mode != "AM":
AMTrackInputTag = cms.InputTag("TTTracksFromPixelDigis", "Level1TTTracks")
	

#from SLHCUpgradeSimulations.L1TrackTrigger.l1TkMuonsExt_cfi import l1TkMuonsExt
#l1TkMuonsExtCSC = l1TkMuonsExt.clone(
#		  L1MuonsInputTag = cms.InputTag("l1extraMuExtended", "csc"),
		  
#	       	  L1TrackInputTag = AMTrackInputTag,
		  #ETAMIN = cms.double(1.1),
#		  )

process.schedule = cms.Schedule(process.p)
#process.l1TkMuonsExtNoZCor = l1TkMuonsExt.clone( correctGMTPropForTkZ = cms.bool(False) )
#process.l1TkMuonsExtCSCNoZCor = l1TkMuonsExtCSC.clone( correctGMTPropForTkZ = cms.bool(False) )
#l1TkMuonsExtSequence = cms.Sequence (l1TkMuonsExt * l1TkMuonsExtCSC
#                                     * l1TkMuonsExtNoZCor * l1TkMuonsExtCSCNoZCor)

#process.schedule = cms.Schedule(process.p, process.TTAssociator_step,process.slhccalo,process.L1Reco,process.pMuons,process.pL1TkPhotons,process.pElectrons,process.pL1TkIsoElectrons,process.pElectronsLoose,process.pAna)
