# Auto generated configuration file
# using: 
# Revision: 1.14 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step2_BE --step DIGI,L1,DIGI2RAW,L1TrackTrigger,RECO:pixeltrackerlocalreco --datatier GEN-SIM-DIGI-RAW --geometry ExtendedPhase2TkBE --beamspot HLLHC --conditions POSTLS261_V2::All --pileup NoPileUp --datamix NODATAMIXER --eventcontent FEVTDEBUG --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1,SLHCUpgradeSimulations/Configuration/phase2TkCustomsBE.customise --filein file:step1P1.root --fileout file:step2BE.root -n 2
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBEReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('file:SingleElectron_Pt5to50_GEN_SIM.root') 
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.14 $'),
    annotation = cms.untracked.string('step2_BE nevts:2'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    fileName = cms.untracked.string('file:SingleElectron_Pt5to50_DIGI_RAW.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW')
    )
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V2::All', '')

#################################################################################################
# Calo trigger information
#################################################################################################
process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff")
process.L1CaloTowerProducer.ECALDigis = cms.InputTag("simEcalTriggerPrimitiveDigis")
process.L1CaloTowerProducer.HCALDigis = cms.InputTag("simHcalTriggerPrimitiveDigis")
process.L1Calo = cms.Path(process.SLHCCaloTrigger)
#################################################################################################
# Beam Spot  
#################################################################################################
process.BeamSpotFromSim=cms.EDProducer("BeamSpotFromSimProducer")
process.BeamSpot  = cms.Path(process.BeamSpotFromSim)
#################################################################################################
# Analizer 
#################################################################################################
process.NtupleMaker = cms.EDAnalyzer('Pxecal')
process.p = cms.Path(process.NtupleMaker)
process.TFileService = cms.Service("TFileService", fileName = cms.string('SingleElectron_Pt5to50_ntuple.root') )

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)
process.reconstruction_step = cms.Path(process.pixeltrackerlocalreco)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.L1TrackTrigger_step,process.reconstruction_step,process.endjob_step,process.FEVTDEBUGoutput_step,process.BeamSpot,process.L1Calo,process.p)

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.phase2TkCustomsBE
from SLHCUpgradeSimulations.Configuration.phase2TkCustomsBE import customise 

#call to customisation function customise imported from SLHCUpgradeSimulations.Configuration.phase2TkCustomsBE
process = customise(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1 

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)

# End of customisation functions
