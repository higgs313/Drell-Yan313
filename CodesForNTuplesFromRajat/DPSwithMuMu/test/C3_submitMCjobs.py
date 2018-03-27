import os,sys
import string, re, shutil
from time import gmtime, localtime, strftime

data = [#"/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",#1
#"/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",#2
#"/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_HCALDebug_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",#3
#"/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM",#4
#"/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-MPI_off/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_ext2-v1/MINIAODSIM", #5
#"/ZZTo4L_13TeV_powheg_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",#5
#"/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",#6
#"/WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",#7
#"/WWToLNuQQ_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",#8
#"/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM",#9
#"/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",#10
#"/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v3/MINIAODSIM",#11
#"/WZJToLLLNu_TuneCUETP8M1_13TeV-amcnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM"#12
"/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/MINIAODSIM"
]

name=[#"TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_18072017",#1
#"WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_18072017",#2
#"DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_18072017",#3
#"DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-herwigpp_18072017",#4
#"DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-MPI_off_18072017",#5
#"ZZTo4L_13TeV_powheg_pythia8_25022017_03032017",#5
#"ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_03032017",#6
#"WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_03032017",#7
#"WWToLNuQQ_13TeV-powheg_03032017",#8
#"WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_03032017",#9
#"WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8_03032017",#10
#"WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_03032017",#11
#"WZJToLLLNu_TuneCUETP8M1_13TeV-amcnlo-pythia8_03032017"#12
"DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp"
]


processname=["HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT","HLT"]
#pileupFile=["S10MC_PUFile.root", "S10MC_PUFile.root", "S10MC_PUFile.root", "S10MC_PUFile.root", "S10MC_PUFile.root", "S10MC_PUFile.root", "S10MC_PUFile.root", "S10MC_PUFile.root", "S7MC_PUFile.root"]

def changeMainConfigFile(outfile,hlt):
    fin  = open("testMC.py")
    pset_cfg      = "py_" + outfile + ".py"
    outfile_root  = "tree_ssWW_13TeV_" + outfile + ".root"
    fout = open(pset_cfg,"w")
    for line in fin.readlines():
        if  line.find("patTuple.root")!=-1:
            line=line.replace("patTuple.root",outfile_root)
        #if  line.find("processname")!=-1:
         #   line=line.replace("processname",processname[hlt]) 
#        if  line.find("pileupFile.root")!=-1:
 #           line=line.replace("pileupFile.root",pileupFile[hlt]) 
        fout.write(line)
    print pset_cfg + " has been written.\n"

def changeCrabTemplateFile(outfile, index):
    fin  = open("crab3TemplateMC.py")
    pset_cfg      = "py_" + outfile + ".py"
    pset_crab     = "crabjob_" + outfile + ".py"
    outfile_root  = "tree_ssWW_13TeV_" + outfile + ".root"
    fout = open(pset_crab,"w")
    for line in fin.readlines():
        if  line.find("mydatapath")!=-1:
            line=line.replace("mydatapath",data[index])
        if  line.find("myanalysis")!=-1:
            line=line.replace("myanalysis",pset_cfg)
        if  line.find("myrootfile")!=-1:
            line=line.replace("myrootfile",outfile_root)
        if  line.find("storePath")!=-1:
            line=line.replace("storePath",outfile)
        if  line.find("workingdir")!=-1:
            line=line.replace("workingdir",outfile)
        fout.write(line)

    print pset_crab + " has been written.\n"

     

###################
import sys

if __name__ == '__main__':

    print sys.argv
if len(sys.argv) <= 1 :
    print "no arguments?"
    print "Usage to submit:     python C3_submitMCjobs.py create"
    print "Usage to get status: python C3_submitMCjobs.py status"
    exit()
if len(sys.argv) > 1 :
    if sys.argv[1] == 'report' :
        for i in range(len(name)):
            reportcommand = "crab report -d " + "crab" + name[i]
            child   = os.system(reportcommand)
    if sys.argv[1] == 'status' :
        for i in range(len(name)):
            statuscommand = "crab status -d " + "crab_" + name[i]
            child   = os.system(statuscommand)
    if sys.argv[1] == 'create' :    
        for i in range(len(name)):
            changeMainConfigFile(name[i],i)
            changeCrabTemplateFile(name[i],i)
        for i in range(len(name)):
            createcommand = "crab submit -c " + "crabjob_" + name[i] + ".py"
            child   = os.system(createcommand)



 #   submitcommand2 = "crab -submit 1-500 -c " + name[i]
  #  child2   = os.system(submitcommand2)
   # submitcommand3 = "crab -submit 501-1000 -c " + name[i]
    #child3   = os.system(submitcommand3)
 #   submitcommand4 = "crab -submit 1001-1500 -c " + name[i]
  #  child4   = os.system(submitcommand4)
   # submitcommand5 = "crab -submit 1501-2000 -c " + name[i]
    #child5   = os.system(submitcommand5)
    #submitcommand6 = "crab -submit 2001-2500 -c " + name[i]
    #child6   = os.system(submitcommand6)
    #submitcommand7 = "crab -submit 2501-3000 -c " + name[i]
    #child7   = os.system(submitcommand7)
    #submitcommand8 = "crab -submit 3001-3500 -c " + name[i]
    #child8   = os.system(submitcommand8)
    #submitcommand9 = "crab -submit 3501-4000 -c " + name[i]
    #child9   = os.system(submitcommand9)
    #submitcommand10 = "crab -submit 4001-4500 -c " + name[i]
    #child10  = os.system(submitcommand10)
    #submitcommand11 = "crab -submit 4501-5000 -c " + name[i]
    #child11  = os.system(submitcommand11)

