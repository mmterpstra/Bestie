version 1.0

import "../structs.wdl" as structs

task CollectMultipleMetrics {
    input {
        File inputBam
        String outputMetricsBasename
        Reference reference
        String gatkModule = "GATK"
        Int? memoryGb = "4"
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 40
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
        ml ~{gatkModule}
        gatk --java-options "-Xmx3g -XX:ParallelGCThreads=1" CollectMultipleMetrics \
        --REFERENCE_SEQUENCE ~{reference.fasta} \
        --INPUT ~{inputBam} \
        --OUTPUT ~{outputMetricsBasename}
    }
    
    output {
        File alignmentMetrics = outputMetricsBasename + "_multiplemetrics" + ".alignment_summary_metrics"
        File baseDistributionMetrics = outputMetricsBasename + "_multiplemetrics" + ".base_distribution_by_cycle_metrics"
        File baseDistributionPdf = outputMetricsBasename + "_multiplemetrics" + ".base_distribution_by_cycle.pdf"
        File insertSizeMetrics = outputMetricsBasename + "_multiplemetrics" + ".insert_size_metrics"
        File insertSizePdf = outputMetricsBasename + "_multiplemetrics" + ".insert_size_histogram.pdf"
        File qualityByCycleMetrics = outputMetricsBasename + "_multiplemetrics" + ".quality_by_cycle_metrics"
        File qualityByCyclePdf = outputMetricsBasename + "_multiplemetrics" + ".quality_by_cycle.pdf"
        File qualityDistributionMetrics = outputMetricsBasename + "_multiplemetrics" + ".quality_distribution_metrics"
        File qualityDistributionPdf = outputMetricsBasename + "_multiplemetrics" + ".quality_distribution.pdf"
        File readLengthPdf = outputMetricsBasename + "_multiplemetrics" + ".read_length_histogram.pdf"
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
    }
}

task CollectHsMetrics {
    input {
        File inputBam
        File targetIntervalList
        String outputMetricsBasename
        Reference reference
        String gatkModule = "GATK"
        Int? memoryGb = "4"
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 40
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
        ml ~{gatkModule}
        gatk --java-options "-Xmx3g -XX:ParallelGCThreads=1" CollectHsMetrics \
        --REFERENCE_SEQUENCE ~{reference.fasta} \
        --BAIT_INTERVALS ~{targetIntervalList} \
        --TARGET_INTERVALS ~{targetIntervalList}\ 
        --INPUT ~{inputBam} \
        --OUTPUT ~{outputMetricsBasename}
    }
    
    output {
        String outHsMetrics = outputMetricsBasename
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
    }
}

#SPLIT_TO_N_READS
#https://github.com/broadinstitute/warp/blob/develop/tasks/broad/Alignment.wdl#L128