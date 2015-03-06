#!/bin/bash


files=`ls /net/scratch_cms/institut_3a/data_Muon_POG/output_histos | grep RecoOnly | grep UT`
for file in $files
do
        python draw.py /net/scratch_cms/institut_3a/data_Muon_POG/output_histos/$file

done
