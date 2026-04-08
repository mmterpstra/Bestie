version 1.0

import "../structs.wdl"
import "../tasks/gatk.wdl" as gatk
import "../tasks/picard.wdl" as picard
import "../tasks/fgbio.wdl" as fgbio
import "../tasks/samtools.wdl" as samtools
import "../tasks/common.wdl" as common

workflow bamQualityControl {
    input {
        File inputBam
        File inputBai
        String gatkModule
        String picardModule
        String fgbioModule
        String samtoolsModule
        String outputPrefix
        Boolean runSamtools=false
        Reference reference
        File? targetIntervalList
        IndexedFile? commonVariants
        Boolean byReadGroup = false
        Boolean flattenArchive = true
    }
    call gatk.CollectMultipleMetrics as collectMultipleMetrics {
        input:
            gatkModule = gatkModule,
            reference = reference,
            inputBam = inputBam,
            outputMetricsBasename = outputPrefix + "_multiplemetrics",
            byReadGroup = byReadGroup
    }
    #this one is sloooow:
    #call gatk.CollectWgsMetrics as wgsMetrics {
    #    input:
    #        gatkModule = gatkModule,
    #        reference = reference,
    #        inputBam = inputBam,
    #        outputMetricsBasename = outputPrefix,
    #}
    if(defined(targetIntervalList)) {
        call gatk.CollectHsMetrics as collectHsMetrics {
            input:
                gatkModule = gatkModule,
                reference = reference,
                targetIntervalList = select_first([targetIntervalList]),
                inputBam = inputBam,
                outputMetricsBasename = outputPrefix,
                byReadGroup = byReadGroup

        }
        call gatk.DepthOfCoverage as depthOfCoverage {
            input:
                gatkModule = gatkModule,
                reference = reference,
                targetIntervalList = select_first([targetIntervalList]),
                inputBam = inputBam,
                outputMetricsBasename = outputPrefix,
        }
        call fgbio.ErrorRateByReadPosition as errorByReadPos {
            input:
                fgbioModule = fgbioModule,
                reference = reference,
                inputBam = inputBam,
                outputBasename = outputPrefix,
                intervals = select_first([targetIntervalList]),
                variants = select_first([commonVariants])
        }     
    }
    if(runSamtools) {
        call samtools.Stats as samStats {
            input:
                samtoolsModule = samtoolsModule,
                inputBam = inputBam,
                inputBai = inputBai,
                reference = reference,
                outputBasename = outputPrefix,
            }
        call samtools.XYIdxStats as samXYIdxStats {
            input:
                samtoolsModule = samtoolsModule,
                inputBam = inputBam,
                inputBai = inputBai,
                outputBasename = outputPrefix,
        }
    }

    call common.ZipFiles as CreateQcZip {
        input:
            fileList = 
                [   collectMultipleMetrics.alignmentMetrics, 
                    collectMultipleMetrics.baseDistributionMetrics,
                    collectMultipleMetrics.baseDistributionPdf,
                    collectMultipleMetrics.insertSizeMetrics,
                    collectMultipleMetrics.insertSizePdf,
                    collectMultipleMetrics.qualityByCycleMetrics,
                    collectMultipleMetrics.qualityByCyclePdf,
                    collectMultipleMetrics.qualityDistributionMetrics,
                    collectMultipleMetrics.readLengthPdf
                    #,wgsMetrics.wgsMetrics
                ],
            optionalFileList = [
                    #dummy file filling
                    select_first([collectHsMetrics.hsMetrics,collectMultipleMetrics.alignmentMetrics]),
                    select_first([depthOfCoverage.dcovMetrics,collectMultipleMetrics.alignmentMetrics]),
                    select_first([errorByReadPos.metrics,collectMultipleMetrics.alignmentMetrics]),
                    select_first([samStats.stats,collectMultipleMetrics.alignmentMetrics]),
                    select_first([samXYIdxStats.idxstats,collectMultipleMetrics.alignmentMetrics])
                ],
            outputPrefix = outputPrefix,
            flattenArchive = flattenArchive
    }
    output {
        #File alignmentMetrics = CollectMultipleMetrics.alignmentMetrics
        #File baseDistributionMetrics = CollectMultipleMetrics.baseDistributionMetrics
        #File baseDistributionPdf = CollectMultipleMetrics.baseDistributionPdf
        #File insertSizeMetrics = CollectMultipleMetrics.insertSizeMetrics
        #File insertSizePdf = CollectMultipleMetrics.insertSizePdf
        #File qualityByCycleMetrics = CollectMultipleMetrics.qualityByCycleMetrics
        #File qualityByCyclePdf = CollectMultipleMetrics.qualityByCyclePdf
        #File qualityDistributionMetrics = CollectMultipleMetrics.qualityDistributionMetrics
        #File qualityDistributionPdf = CollectMultipleMetrics.qualityDistributionPdf
        #File readLengthPdf = CollectMultipleMetrics.readLengthPdf
        #File? hsMetrics = select_first([CollectHsMetrics.hsMetrics,outputPrefix + ".hs_metrics_skipped"])
        File qcZip = CreateQcZip.zip
    }
}