version 1.0

import "../structs.wdl"
import "../tasks/gatk.wdl" as gatk
import "../tasks/picard.wdl" as picard
import "../tasks/common.wdl" as common

workflow bamQualityControl {
    input {
        File inputBam
        File inputBai
        String gatkModule
        String picardModule
        String outputPrefix
        Reference reference
        File? targetIntervalList
    }
    call gatk.CollectMultipleMetrics as CollectMultipleMetrics {
        input:
            gatkModule = gatkModule,
            reference = reference,
            inputBam = inputBam,
            outputMetricsBasename = outputPrefix + "_multiplemetrics"
    }
    if(defined(targetIntervalList)) {
        call gatk.CollectHsMetrics as CollectHsMetrics {
            input:
                gatkModule = gatkModule,
                reference = reference,
                targetIntervalList = select_first([targetIntervalList]),
                inputBam = inputBam,
                outputMetricsBasename = outputPrefix
        }
    }
    call common.ZipFiles as CreateQcZip {
        input:
            fileList = 
                [   CollectMultipleMetrics.alignmentMetrics, 
                    CollectMultipleMetrics.baseDistributionMetrics,
                    CollectMultipleMetrics.baseDistributionPdf,
                    CollectMultipleMetrics.insertSizeMetrics,
                    CollectMultipleMetrics.insertSizePdf,
                    CollectMultipleMetrics.qualityByCycleMetrics,
                    CollectMultipleMetrics.qualityByCyclePdf,
                    CollectMultipleMetrics.qualityDistributionMetrics,
                    CollectMultipleMetrics.readLengthPdf
                ],
            optionalFileList = [select_first([CollectHsMetrics.hsMetrics,outputPrefix + ".hs_metrics_skipped"])],
            outputPrefix = outputPrefix + "_qc"
    }
    output {
        File alignmentMetrics = CollectMultipleMetrics.alignmentMetrics
        File baseDistributionMetrics = CollectMultipleMetrics.baseDistributionMetrics
        File baseDistributionPdf = CollectMultipleMetrics.baseDistributionPdf
        File insertSizeMetrics = CollectMultipleMetrics.insertSizeMetrics
        File insertSizePdf = CollectMultipleMetrics.insertSizePdf
        File qualityByCycleMetrics = CollectMultipleMetrics.qualityByCycleMetrics
        File qualityByCyclePdf = CollectMultipleMetrics.qualityByCyclePdf
        File qualityDistributionMetrics = CollectMultipleMetrics.qualityDistributionMetrics
        File qualityDistributionPdf = CollectMultipleMetrics.qualityDistributionPdf
        File readLengthPdf = CollectMultipleMetrics.readLengthPdf
        File? hsMetrics = select_first([CollectHsMetrics.hsMetrics,outputPrefix + ".hs_metrics_skipped"])
        File qcZip = CreateQcZip.zip
    }
}