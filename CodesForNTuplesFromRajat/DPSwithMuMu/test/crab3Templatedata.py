from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'workingdir'
config.General.workArea = '/afs/cern.ch/work/r/rgupta/CMSSW_8_0_20/src/ZplusJets/DPSwithMuMu/test'
config.General.transferOutputs = True
config.General.transferLogs = True


config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'myanalysis'
#config.JobType.inputFiles = ['MuScleFit_2012ABC_DATA_ReReco_53X.txt', 'MuScleFit_2012D_DATA_53X.txt', 'MuScleFit_2012D_DATA_ReReco_53X.txt', 'MuScleFit_2012_MC_53X_smearReReco.txt']
#config.JobType.outputFiles = ['myrootfile']

config.Data.inputDataset = 'mydatapath'

config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 150000
config.Data.totalUnits = -1 #number of event
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'

config.Data.outLFNDirBase = '/store/user/rgupta/DPSinZpJ/04032017/'
config.Data.publication = False
config.Data.outputDatasetTag = 'workingdir'
config.Site.storageSite = 'T2_IN_TIFR'

