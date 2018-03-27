from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'workingdir'
#config.General.workArea = '/afs/cern.ch/work/a/anmehta/work/SL6_setup/CMSSW_5_3_20/src/TestPAT/MPIwithWW/CRAB3'
config.General.workArea = '/afs/cern.ch/work/r/rgupta/CMSSW_8_0_20/src/ZplusJets/DPSwithMuMu/test'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'myanalysis'
#config.JobType.inputFiles = ['MuScleFit_2012ABC_DATA_ReReco_53X.txt', 'MuScleFit_2012D_DATA_53X.txt', 'MuScleFit_2012D_DATA_ReReco_53X.txt', 'MuScleFit_2012_MC_53X_smearReReco.txt']
#config.JobType.outputFiles = ['myrootfile']


config.Data.inputDataset = 'mydatapath'
config.Data.splitting    = 'FileBased'  #'LumiBased'
config.Data.unitsPerJob  = 3 

config.Data.outLFNDirBase = '/store/user/rgupta/DPSinZpJ/18072017/'
config.Data.publication = False
config.Data.outputDatasetTag = 'workingdir'
config.Site.storageSite = 'T2_IN_TIFR'

#config.section_('User')
#config.section_('Site')

