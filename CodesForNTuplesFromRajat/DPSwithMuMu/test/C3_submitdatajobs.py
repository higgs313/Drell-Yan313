import os,sys
import string, re, shutil
from time import gmtime, localtime, strftime
data = ["/SingleMuon/Run2016B-23Sep2016-v3/MINIAOD",
	"/SingleMuon/Run2016C-23Sep2016-v1/MINIAOD",
	"/SingleMuon/Run2016D-23Sep2016-v1/MINIAOD",
	"/SingleMuon/Run2016E-23Sep2016-v1/MINIAOD",
	"/SingleMuon/Run2016F-23Sep2016-v1/MINIAOD",
	"/SingleMuon/Run2016G-23Sep2016-v1/MINIAOD",
	"/SingleMuon/Run2016H-PromptReco-v3/MINIAOD",
	"/SingleMuon/Run2016H-PromptReco-v2/MINIAOD",
	#"/SingleMuon/Run2016H-PromptReco-v1/MINIAOD"
        ]

name=["SingleMuon_Run2016B-23Sep2016-v3_18072017_MINIAOD",
      "SingleMuon_Run2016C-23Sep2016-v1_18072017_MINIAOD",
      "SingleMuon_Run2016D-23Sep2016-v1_18072017_MINIAOD",
      "SingleMuon_Run2016E-23Sep2016-v1_18072017_MINIAOD",
      "SingleMuon_Run2016F-23Sep2016-v1_18072017_MINIAOD",
      "SingleMuon_Run2016G-23Sep2016-v1_18072017_MINIAOD",
      "SingleMuon_Run2016H-PromptReco-v3_18072017_MINIAOD",
      "SingleMuon_Run2016H-PromptReco-v2_18072017_MINIAOD",
      #"SingleMuon_Run2016H-PromptReco-v1_17012017_MINIAOD"
    ]




processname=["HLT","HLT","HLT","HLT"]


def changeMainConfigFile(outfile,hlt):
    fin  = open("testData.py")
    pset_cfg      = "py_" + outfile + ".py"
    outfile_root  = "tree_ssWW_13TeV_" + outfile + ".root"
    fout = open(pset_cfg,"w")
    for line in fin.readlines():
        if  line.find("patTuple.root")!=-1:
            line=line.replace("patTuple.root",outfile_root)
        fout.write(line)
    print pset_cfg + " has been written.\n"

def changeCrabTemplateFile(outfile, index):
    fin  = open("crab3Templatedata.py")
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
for i in range(len(name)):
    changeMainConfigFile(name[i],i)
    changeCrabTemplateFile(name[i],i)

for i in range(len(name)):
    createcommand = "crab submit -c " + "crabjob_" + name[i] + ".py"
    child   = os.system(createcommand)


