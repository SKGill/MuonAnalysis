from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
config.General.requestName   = 'Run2017_MC_p100'
config.General.transferLogs = True
config.General.workArea = 'crab'
config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName    = 'ntuple.py'
#config.JobType.inputFiles = ['startup-v1_DESRUN2_74_V4_ape-candidate1.db','startup-v1_DESRUN2_74_V4_ape-candidate2.db']

config.section_("Data")
config.Data.inputDataset = '/SPLooseMuCosmic_38T_p500/RunIISummer17CosmicDR-94X_mc2017cosmics_realistic_deco_v3-v1/GEN-SIM-RECO'
# '/SPLooseMuCosmic_38T_p100-500/RunIISummer17CosmicDR-94X_mc2017cosmics_realistic_deco_v3-v1/GEN-SIM-RECO'
#'/SPLooseMuCosmic_38T_p500/RunIISummer17CosmicDR-94X_mc2017cosmics_realistic_deco_v3-v1/GEN-SIM-RECO'
#config.Data.allowNonValidInputDataset = True
#f = open('cosmic_files.txt')
#flist = []
#for line in f:
 #   flist.append(line)
#config.Data.userInputFiles = flist
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.publication = True
# This string is used to construct the output dataset name
#config.Data.publishDataName = 'CRAB3-tutorial'
config.Data.outLFNDirBase = '/store/group/phys_smp/skaur/MC_2017'
#config.Data.ignoreLocality = True

# These values only make sense for processing data
#    Select input data based on a lumi mask
#config.Data.lumiMask = 'cosmics15_CRAFT_238361_239517_json.txt'
#    Select input data based on run-ranges
#config.Data.runRange = '190456-194076'

config.section_("Site")
# Where the output files will be transmitted to
config.Site.storageSite = 'T2_CH_CERN'

#config.Site.whitelist = ['T2_US*']
