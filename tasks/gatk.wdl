version 1.0

import "../structs.wdl" as structs

task CollectMultipleMetrics {
    input {
        File inputBam
        String outputMetricsBasename
        Reference reference
        String gatkModule = "GATK"
        Int? memoryGb = "4"
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
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
        File alignmentMetrics = outputMetricsBasename + ".alignment_summary_metrics"
        File baseDistributionMetrics = outputMetricsBasename + ".base_distribution_by_cycle_metrics"
        File baseDistributionPdf = outputMetricsBasename + ".base_distribution_by_cycle.pdf"
        File insertSizeMetrics = outputMetricsBasename + ".insert_size_metrics"
        File insertSizePdf = outputMetricsBasename + ".insert_size_histogram.pdf"
        File qualityByCycleMetrics = outputMetricsBasename + ".quality_by_cycle_metrics"
        File qualityByCyclePdf = outputMetricsBasename + ".quality_by_cycle.pdf"
        File qualityDistributionMetrics = outputMetricsBasename + ".quality_distribution_metrics"
        File qualityDistributionPdf = outputMetricsBasename + ".quality_distribution.pdf"
        File readLengthPdf = outputMetricsBasename + ".read_length_histogram.pdf"
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
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
        ml ~{gatkModule}
        gatk --java-options "-Xmx3g -XX:ParallelGCThreads=1" CollectHsMetrics \
        --REFERENCE_SEQUENCE "~{reference.fasta}" \
        --BAIT_INTERVALS "~{targetIntervalList}" \
        --TARGET_INTERVALS "~{targetIntervalList}" \
        --INPUT "~{inputBam}" \
        --OUTPUT "~{outputMetricsBasename}.hs_metrics"
    }
    
    output {
        File hsMetrics = outputMetricsBasename + ".hs_metrics"
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
    }
}

#wip
task BaseQualityScoreRecalibration {
    input {
        File inputBam
        File targetIntervalList
        String outputMetricsBasename
        Reference reference
        String gatkModule = "GATK"
        Int? memoryGb = "4"
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
        ml ~{gatkModule}
        gatk --java-options "-Xmx3g -XX:ParallelGCThreads=1" BaseRecalibrator \
        --REFERENCE_SEQUENCE "~{reference.fasta}" \
        --BAIT_INTERVALS "~{targetIntervalList}" \
        --TARGET_INTERVALS "~{targetIntervalList}" \
        --INPUT "~{inputBam}" \
        --OUTPUT "~{outputMetricsBasename}.hs_metrics"
    }
    
    output {
        File hsMetrics = outputMetricsBasename + ".hs_metrics"
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
    }
}

task HaplotypeCaller {
    input {
        File inputBam
        File inputBai
        File targetIntervalList
        String outputVcfBasename
        Reference reference
        String gatkModule = "GATK"
        Int memoryGb = "4"
        Int javaMemoryGb = memoryGb - 1
        Boolean makeGvcf = true
        Boolean makeBamOut = true
        Boolean useSpanningEventGenotyping = false
        Float? contamination = 0
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
    }

    String vcfSuffix =  if makeGvcf then ".g.vcf" else ".vcf"
    String bamoutArg =  if makeBamOut then "-bamout " + outputVcfBasename + ".bamout.bam" else ""
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
        set -e

        ml ~{gatkModule}
        gatk --java-options "-Xmx${javaMemoryGb}g -Xms${javaMemoryGb}g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        HaplotypeCaller \
        -R ~{reference.fasta} \
        -I ~{inputBam} \
        -L ~{targetIntervalList} \
        -O "~{outputVcfBasename}~{vcfSuffix}" \
        -contamination ~{default=0 contamination} \
        -G StandardAnnotation -G StandardHCAnnotation ~{true="-G AS_StandardAnnotation" false="" makeGvcf} \
        ~{false="--disable-spanning-event-genotyping" true="" useSpanningEventGenotyping} \
        -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
        ~{true="-ERC GVCF" false="" makeGvcf} \
        ~{bamoutArg}

        touch ~{outputVcfBasename}.bamout.bam
    }
    
    output {
        File outVcf = outputVcfBasename + vcfSuffix
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
    }
}

#SPLIT_TO_N_READS
#https://github.com/broadinstitute/warp/blob/develop/tasks/broad/Alignment.wdl#L128