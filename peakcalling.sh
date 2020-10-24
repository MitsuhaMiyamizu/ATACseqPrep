#!/bin/bash
#Peak Calling

ls .bed | while read id ;do (macs2 callpeak -t $id --call-summits -B -q 0.05 -g hs --nomodel -n ${id%%.} --outdir ../peaks/) ;done
macs2 callpeak -t example.filtered.bed -q 0.05 -g hs --nomodel --call-summits -B --extsize 200 -n example --outdir ../peaks/
