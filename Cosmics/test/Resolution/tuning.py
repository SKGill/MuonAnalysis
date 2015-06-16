#!/usr/bin/env python

import os,sys
#sys.argv.append('-b')
import ROOT

# 0,1 Global
# 2,3 StandAlone
# 4,5 TkOnly
# 6,7 TPFMS
# 8,9 Picky
# 10,11 DYT
a={'0':2,'1':2,'4':1,'5':1,'6':3,'7':3,'8':4,'9':4,'10':5,'11':5}
b={'0':1,'1':1,'4':2,'5':2,'6':3,'7':3,'8':4,'9':4,'10':5,'11':5}


def Dptcut(event,cut,j):
    dptmin = 1
    if cut > 0:
        for i in [0,2,3,4,5]:
            if event.unprop_error_pt[2*i+j]/event.unprop_pt[2*i+j] < 1:
                dptmin = event.unprop_error_pt[2*i+j]/event.unprop_pt[2*i+j]
    return dptmin
        


ROOT.gStyle.SetOptStat(111111)
f = ROOT.TFile('~/nobackup/resolution_ntuple_v54_newalign.root')
d = f.Get('UTpickedTracks')
t = d.Get('t')
run_list = [239213, 239286, 239513, 238747, 239283, 239514, 239515, 238983, 238981, 239044, 239203, 239200, 238985, 238989, 238474, 238486, 238643, 238642, 239179, 238649, 239509, 239425, 239291, 239359, 239507, 239517, 239215, 239214, 239194, 239195, 239192, 238746, 239058, 238749, 239219, 239199, 239207, 239160, 238671, 239342, 239340, 239368, 239201, 238592, 238515, 238514, 239186, 238511, 238513, 238512, 238757, 239020, 238519, 238752, 238751, 238750, 239362, 238668, 238663, 238662, 238660, 238843, 238443, 239404, 238445, 238685, 238522, 238723, 238720, 238721, 238724, 239230, 239334, 238832, 238830, 238831, 238835, 239004, 239264, 239413, 238637, 238864, 238828, 238825, 239485, 238826, 239399, 238708, 238709, 238542, 238547, 238545, 238696, 238810, 239184, 238818, 239224, 238713, 238716, 239262, 238491, 238361, 238495, 239060, 239061, 239305, 238556, 239387, 238807, 239386, 239189, 239150, 238875, 238780, 239278, 239221, 238764, 239371, 238876, 238877, 238799, 238505, 239382, 239211, 238972]

# Run tuneP from ntuple, ie, rewrite the cocktail
h = ROOT.TH1F('h','',5,1,6)
g = ROOT.TH1F('g','',100,-0.1,4000)

for event in t:
    if event.run not in run_list:
        continue
    pairs = [8,9] # Start with picky
    for j in [0,1]:
        dptcut = 0.25
        dptmin = Dptcut(event,dptcut,j)
        if dptmin > dptcut:
            dptcut = dptmin+0.15
        probs = [0,0,0,0,0]        
        for x,i in enumerate([2,0,3,4,5]): # TkOnly,Global,TPFMS,Picky,DYT
            if (event.pixel_hits[2*i+j] + event.strip_hits[2*i+j]) > 0 and (event.unprop_error_pt[2*i+j]/event.unprop_pt[2*i+j] < dptcut or dptcut < 0):
                probs[x] = event.prob[2*i+j]

        if probs[3] == 0:
            if dptcut > 0:
                if probs[4] > 0.:
                    pairs[j] = 10+j
                elif probs[0] > 0.:
                    pairs[j] = 4+j
                elif probs[2] > 0.:
                    pairs[j] = 6+j
                elif probs[1] > 0.:
                    pairs[j] = j
                    
        if probs[4] > 0 and probs[3] > 0:            
            if event.unprop_pt[8+j] > 0 and event.unprop_pt[10+j] > 0 and (event.unprop_error_pt[10+j]/event.unprop_pt[10+j] - event.unprop_error_pt[8+j]/event.unprop_pt[8+j]) <= 0:
                pairs[j] = 10+j        
        if probs[0] > 0 and probs[a[str(pairs[j])]-1] > 0 and (probs[a[str(pairs[j])]-1] - probs[0]) > 17:
            pairs[j] = 4+j
        if probs[2] > 0 and (probs[a[str(pairs[j])]-1] - probs[2]) > 40:
            pairs[j] = 6+j
        if event.unprop_pt[pairs[j]] < 200 or event.unprop_pt[4+j] < 200:
            pairs[j] = 4+j                            

    h.Fill(b[str(pairs[0])])
    h.Fill(b[str(pairs[1])])

    qinvpt_0 = event.charge[pairs[0]]*1./event.pt[pairs[0]]
    qinvpt_1 = event.charge[pairs[1]]*1./event.pt[pairs[1]]
    #g.Fill((qinvpt_0-qinvpt_1)/qinvpt_1)
    #print event.event,event.unprop_pt[pairs[1]]-event.pt[pairs[1]]
    g.Fill((event.pt[pairs[1]]))
    
h.Draw()

            
