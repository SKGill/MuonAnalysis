from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
config.General.requestName   = 'Muon_Resolution_Ntuple_FR_R_74_V15B_tuneR_conservative'
config.General.transferLogs = True
config.General.workArea = 'crab'
config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName    = 'ntuple.py'
config.JobType.inputFiles = ['startup-v1_DESRUN2_74_V4_ape-candidate1.db','startup-v1_DESRUN2_74_V4_ape-candidate2.db']

config.section_("Data")
#config.Data.inputDataset = '/Cosmics/Commissioning2015-CosmicSP-04Jun2015-v1/RAW-RECO'
#config.Data.allowNonValidInputDataset = True
f = open('cosmic_files.txt')
flist = []
for line in f:
    flist.append(line)
config.Data.userInputFiles = flist
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.publication = True
# This string is used to construct the output dataset name
#config.Data.publishDataName = 'CRAB3-tutorial'
config.Data.outLFNDirBase = '/store/user/jchavesb/Muon_Resolution_Ntuple_FR_R_74_V15B_tuneR_conservative'
config.Data.ignoreLocality = True

# These values only make sense for processing data
#    Select input data based on a lumi mask
#config.Data.lumiMask = 'cosmics15_CRAFT_238361_239517_json.txt'
#    Select input data based on run-ranges
#config.Data.runRange = '190456-194076'

config.section_("Site")
# Where the output files will be transmitted to
config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.whitelist = ['T2_US*']
