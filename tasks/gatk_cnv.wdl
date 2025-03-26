version 1.0

import "../structs.wdl" as structs

#tbd
#PreprocessIntervals -> (run for tumor/normal)[ CollectCounts -> CollectAllelicCounts -> DenoiseReadCounts  \
#-> ModelSegments -> CallCopyRatioSegments -> PlotDenoisedCopyRatios -> PlotModeledSegments ] -> Annotatewith funco/oncotator
#resource https://github.com/gatk-workflows/gatk4-somatic-cnvs/blob/master/cnv_somatic_pair_workflow.wdl
