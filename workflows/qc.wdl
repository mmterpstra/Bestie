version 1.0

import "../structs.wdl"
import "../tasks/gatk.wdl" as gatk
import "../tasks/picard.wdl" as picard

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
                outputMetricsBasename = outputPrefix + "_hsmetrics"
        }
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
    }
}