#!/usr/bin/env python

import argparse, os, pprint, itertools, FWCore.ParameterSet.Config as cms
from FWCore.PythonUtilities.LumiList import LumiList
from bins import make_bins
from runs import get_run_list
########################################################################################

# Parse options.

# Our options are prefixed with + instead of - or -- because cmsRun eats
# options starting with dashes, and even throws an error when it finds
# one it doesn't recognize.
parser = argparse.ArgumentParser(prefix_chars='+', description='CosmicSplittingResolutionHistos (histogram maker)')

# The filename is the required positional argument. This will also
# consume cmsRun arguments like -j and the python script filename.
# (argparse will treat them as positionals since they don't start with
# the prefix char '+'.)
parser.add_argument('input', nargs='+',
                    help='The filename for the input ROOT file containing the ntuple.')

parser.add_argument('+input-directory', default='UTpickedTracks',
                    help='In the input ROOT file, the directory the ntuple TTree "t" is in. Default is "%(default)s".')
parser.add_argument('+output-fn',
                    help='Override the output filename. (Default is to write to the current directory a file with basename the same as the input, but ending in .histos.root.)')
parser.add_argument('+output-fn-tag',
                    help='A tag for the output filename.')
parser.add_argument('+no-use-run-list', action='store_false', dest='use_run_list',
                    help='For data, ignore any run selection criteria and use all runs in the input file.')
parser.add_argument('+run-type', action='append', dest='run_types', choices='cosmics collisions commissioning'.split(),
                    help='For data, include runs of type TYPE. (May be specified more than once; default is %(default)s.)')
parser.add_argument('+require-rpc-good', action='store_true',
                    help='In the run selection, require RPC subdetector marked as good in DQM.')
parser.add_argument('+min-run', type=int, default=0,
                    help='For data, drop all runs below MIN_RUN.')
parser.add_argument('+max-run', type=int, default=9999999,
                    help='For data, drop all runs below MIN_RUN.')
parser.add_argument('+is-mc', action='store_true',
                    help='Specified input file is MC. (Data assumed by default.)')
parser.add_argument('+only-sample', type=int, default=-1,
                    help='Only use events from the specified sample id (0 = data, 1-3 = MC for p > 10, p > 100, p > 500 samples). -1 = default = use anything.')
parser.add_argument('+no-wrong-sample', action='store_false', dest='check_for_wrong_sample',
                    help='Do not check for "wrong" sample, i.e. use all MC events regardless of dataset id.')
parser.add_argument('+copy-selected-events', action='store_true',
                    help='Copy events passing the overall selection (in the good run list, numbers of track hits, etc.) to a new ntuple in the output file.')
parser.add_argument('+min-pixel-layers', type=int, default=1,
                    help='Require MIN_PIXEL_LAYERS pixel layers on every track (except standalone). Default is %(default)s.')
parser.add_argument('+min-strip-layers', type=int, default=5,
                    help='Require MIN_STRIP_LAYERS strip layers on every track (except standalone). Default is %(default)s.')
parser.add_argument('+max-dxy', type=float, default=100,
                    help='')
parser.add_argument('+max-dz', type=float, default=100,
                    help='')
parser.add_argument('+min-muon-hits', type=int, default=1,
                    help='Require MIN_MUON_HITS muon hits on every track (except tracker-only). Default is %(default)s.')
parser.add_argument('+no-check-shared-hits', action='store_false', dest='check_shared_hits',
                    help='Do not check whether there are shared hits between upper, lower tracks.')
parser.add_argument('+muon-subdet', default='dtonly', choices='dtonly either'.split(),
                    help='Which muon subdetectors to allow (and require DQM flags good for). Choices are %(choices)s; default is "%(default)s".')
parser.add_argument('+require-tt25', action='store_true',
                    help='Require technical trigger bit 25 to have fired in each event.')
parser.add_argument('+require-not-tt25', action='store_true',
                    help='Require technical trigger bit 25 not to have fired in each event.')
parser.add_argument('+bin-by-run', action='store_true',
                    help='Instead of binning by pT (default), bin by run number.')
parser.add_argument('+bin-by-phi', action='store_true',
                    help='Instead of binning by pT (default), bin by phi.')
parser.add_argument('+make-scatterplots', action='store_true',
                    help='Also make the huge TH2 scatterplots.')
options = parser.parse_args()
#JMTBAD groups

################################################################################

# Finalize the options after including any overrides, and do some
# basic checks of consistency.

options.input_fn = [x for x in options.input if x.endswith('.root')]
if len(options.input_fn) != 1:
    raise argparse.ArgumentError('one positional argument required: input_filename.root')
options.input_fn = options.input_fn[0]

if not options.output_fn:
    options.output_fn = os.path.basename(options.input_fn).replace('.root', '.histos.root')

if options.output_fn_tag:
    options.output_fn = options.output_fn.replace('.histos', '_' + options.output_fn_tag + '.histos')
    
if not options.run_types:
    options.run_types = ['cosmics']

options.require_pixels = options.min_pixel_layers > 0

if options.require_tt25 and options.require_not_tt25:
    raise argparse.ArgumentError('cannot require_tt25 and require_not_tt25')
if options.require_tt25 and not options.require_rpc_good:
    print 'warning: require_tt25 and not require_rpc_good!'

################################################################################

run_list = []
if options.use_run_list and not options.is_mc:
    run_list = []
    #run_list = [r for r in get_run_list(options.run_types, options.muon_subdet, options.require_pixels, options.require_rpc_good) if options.min_run <= r <= options.max_run]
    #run_list = [239213, 239286, 239513, 238747, 239283, 239514, 239515, 238983, 238981, 239044, 239203, 239200, 238985, 238989, 238474, 238486, 238643, 238642, 239179, 238649, 239509, 239425, 239291, 239359, 239507, 239517, 239215, 239214, 239194, 239195, 239192, 238746, 239058, 238749, 239219, 239199, 239207, 239160, 238671, 239342, 239340, 239368, 239201, 238592, 238515, 238514, 239186, 238511, 238513, 238512, 238757, 239020, 238519, 238752, 238751, 238750, 239362, 238668, 238663, 238662, 238660, 238843, 238443, 239404, 238445, 238685, 238522, 238723, 238720, 238721, 238724, 239230, 239334, 238832, 238830, 238831, 238835, 239004, 239264, 239413, 238637, 238864, 238828, 238825, 239485, 238826, 239399, 238708, 238709, 238542, 238547, 238545, 238696, 238810, 239184, 238818, 239224, 238713, 238716, 239262, 238491, 238361, 238495, 239060, 239061, 239305, 238556, 239387, 238807, 239386, 239189, 239150, 238875, 238780, 239278, 239221, 238764, 239371, 238876, 238877, 238799, 238505, 239382, 239211, 238972]
    #run_list = [233228, 233235, 233238, 233907, 233947, 233950, 233953, 233957, 233960, 233964, 233967, 233968, 233985, 233989, 233990, 233997, 233999, 234000, 234002, 234020, 234026, 234029, 234033, 234034, 234035, 234037, 234041, 234047, 234052, 234192, 234193, 234195, 234199, 234204, 234205, 234206, 234224, 234230, 234232, 234237, 234239, 234278, 234279, 234284, 234292, 234294, 234297, 234304, 234332, 234343, 234366, 234374, 234389, 234390, 234393, 234397, 234402, 234403, 234404, 234405, 234407, 234410, 234413, 234498, 234500, 234510, 234518, 234530, 234533, 234535, 234541, 234542, 234555, 234556, 234561, 234623, 234628, 234632, 234633, 234634, 234641, 234645, 234646, 234647, 234649, 234650, 234651, 234846, 234848, 234856, 234858, 234873, 234874, 234882, 235050, 235058, 235060, 235062, 235065, 235068, 235073, 235074, 235109, 235130, 235137, 235167, 235172, 235175, 235177, 235178, 235179, 235183, 235185, 235186, 235191, 235192, 235214, 235222, 235241, 235243, 235251, 235253, 235254, 235272, 235278, 235284, 235286, 235295, 235313, 235316, 235328, 235329, 235330, 235346, 235490, 235494, 235495, 235510, 235519, 235523, 235526, 235528, 235686, 235691, 235698, 235699, 235700, 235702, 235707, 235802, 235837, 235840, 235843, 235848, 235849, 235990, 235992, 235994, 235995, 235996, 235997, 235998, 236002, 236154, 236155, 236157, 236159, 236160, 236161, 236164, 236165, 236183, 236185, 236187, 236190, 236191, 236194, 236201, 236203, 236206, 236208, 236209, 236211, 236231, 236232, 236234, 236241, 236250, 236281, 236310, 236312, 236315, 236339, 236351, 236355, 236356, 236487, 236490, 236497, 236516, 236524, 236541, 236542, 236545, 236549, 236551, 236743, 236755, 236763, 236764, 236765, 236767, 236768, 236770, 236886, 236890, 236899, 236905, 237068, 237143, 237146, 237201, 237202, 237203, 237205, 237215, 237216, 237228, 237229, 237230, 237232, 237234, 237239, 237243, 237244, 237245, 237246, 237251, 237255, 237256, 237257, 237260, 237262, 237265, 237315, 237318, 237331, 237332, 237343, 237346, 237348, 237426, 237430, 237431, 237451, 237464, 237465, 237481, 237483, 237485, 237490, 237493, 237588, 237592, 237595, 237602, 237605, 237611, 237614, 237615, 237667, 237700, 237711, 237712, 237715, 237723, 237736, 237738, 237744, 237745, 237748, 237750, 237756, 237761, 237764, 237869, 237890, 237893, 237895, 237897, 237898, 237901, 237904, 237909, 237913, 237916, 237919, 237922, 237925, 237929, 237945, 237951, 237952, 237953, 237954, 237955, 237956, 237958, 237962, 237963, 237967, 238076, 238077, 238081, 238194, 238198, 238199, 238201, 238203, 238205, 238206, 238207, 238208, 238209, 238210, 238211, 238283, 238285, 238309, 238317, 238318, 238329, 238335, 238484, 238492, 238507, 238517, 238521, 238529, 238711, 238767, 238769, 238771, 238774, 238775, 238862, 238984, 239216, 239284, 239338, 239483]    

l = [x for x in LumiList('cosmics_Run2015_all_pix_strip_DT_RPC_complete.json').getLumis()]
lumi_list = []
#for x in l:
#    lumi_list.append(x[0])
#    lumi_list.append(x[1])

f_ev = open('same_ev.txt')
ev_list = []
#for l in f_ev:
    #ev_list = [int(x) for x in l.strip().split(',')]

cfg = cms.PSet(
    directory               = cms.string(options.input_directory),
    is_mc                   = cms.bool(options.is_mc),
    only_sample             = cms.int32(options.only_sample),
    check_for_wrong_sample  = cms.bool(options.check_for_wrong_sample),
    filename                = cms.string(options.input_fn),
    min_muon_hits           = cms.int32(options.min_muon_hits),
    min_pixel_layers        = cms.int32(options.min_pixel_layers),
    min_strip_layers        = cms.int32(options.min_strip_layers),
    max_tpfms_station       = cms.int32(-1),
    max_dxy                 = cms.double(options.max_dxy),
    max_dz                  = cms.double(options.max_dz),
    no_dt_allowed           = cms.bool(options.muon_subdet == 'csconly'),
    no_csc_allowed          = cms.bool(options.muon_subdet == 'dtonly'),
    check_tksta_dphi        = cms.bool(False),
    check_shared_hits       = cms.bool(options.check_shared_hits),
    use_unpropagated_values = cms.bool(True),
    pp_reco_mode            = cms.bool(False),
    force_run_list          = cms.vuint32(run_list),
    force_lumi_list         = cms.vuint32(lumi_list),
    force_event_list        = cms.vuint32(ev_list),
    require_tt25            = cms.bool(options.require_tt25),
    require_not_tt25        = cms.bool(options.require_not_tt25),
    min_bfield              = cms.double(3.7),
    copy_selected_events    = cms.bool(options.copy_selected_events),
    make_scatterplots       = cms.bool(options.make_scatterplots),
    )

if options.bin_by_run:
    if not run_list:
        raise ValueError('bin_by_run specified but run_list is empty')
    bins = make_bins('run', run_list)
elif options.bin_by_phi:
    bins = make_bins('phi', [x/100. for x in range(-320, 321, 32)])
else:
    bins = make_bins('pt')

print 'configuring config:'

z = cfg.parameters_().items() + [
    ('bins, number of', len(bins)),
    ('bins', bins),
    ('run_types', options.run_types),
    ('require_rpc_good', options.require_rpc_good),
    ('output_file', options.output_fn),
    ]
z.sort(key=lambda x: x[0])
for k,v in z:
    print '%s: %s' % (k.ljust(25), pprint.pformat(v))

########################################################################################

process = cms.Process('CosmicSplittingResolutionHistos')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))
process.source = cms.Source('EmptySource')
process.TFileService = cms.Service('TFileService', fileName=cms.string(options.output_fn))
process.histos = cms.EDAnalyzer('CosmicSplittingResolutionHistos', cfg, bins = cms.VPSet(*[b.pset() for b in bins]))
process.p = cms.Path(process.histos)
