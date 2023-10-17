version 1.0

import "../structs.wdl" as structs

task CollectMultipleMetrics {
    input {
        File inputBam
        String outputMetricsBasename
        Reference reference
        String gatkModule = "GATK"
        Int? memoryGb = "4"
        Int? javaXmxMemoryMb = floor(memoryGb*0.9*1024)
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
        Int disk = ceil(size(inputBam, "M")*1.1)
        Boolean byReadGroup = false
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
        ml ~{gatkModule}
        gatk --java-options "-Xmx~{javaXmxMemoryMb}m -XX:ParallelGCThreads=1" CollectMultipleMetrics \
        --REFERENCE_SEQUENCE ~{reference.fasta} \
        --INPUT ~{inputBam} \
        --OUTPUT ~{outputMetricsBasename} \
        ~{if byReadGroup then "--METRIC_ACCUMULATION_LEVEL READ_GROUP " else ""} 

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
        Int? javaXmxMemoryMb = floor(memoryGb*0.9*1024)
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
        Int disk = ceil(size(inputBam, "M")*2.1)
        Boolean byReadGroup = false
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
        ml ~{gatkModule}

        gatk --java-options "-Xmx~{javaXmxMemoryMb}m -XX:ParallelGCThreads=1" CollectHsMetrics \
        --REFERENCE_SEQUENCE "~{reference.fasta}" \
        --BAIT_INTERVALS "~{targetIntervalList}" \
        --TARGET_INTERVALS "~{targetIntervalList}" \
        --INPUT "~{inputBam}" \
        --OUTPUT "~{outputMetricsBasename}.hs_metrics" \
        ~{if byReadGroup then "--METRIC_ACCUMULATION_LEVEL READ_GROUP " else ""} 

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
        String outputRecalibrationReport
        String gatkModule = "GATK"
        #File targetIntervalList
        #String outputMetricsBasename
        Reference reference
        IndexedFile dbsnp
        Array [IndexedFile] knownSites
        String gatkModule = "GATK"
        Int? memoryGb = "6"
        Int? javaXmxMemoryMb = floor(memoryGb*0.9*1024)
        Int disk = ceil(size(inputBam, "M")*1.1)
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
        #rray [File] knownSitesVcfs = select_all(knownSites).file
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command <<<
        ml ~{gatkModule}
        KNOWNSITES=$(cat ~{write_objects(knownSites)}| \
            cut -f 1 | \
            tail -n+2)
        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms5000m -Xmx~{javaXmxMemoryMb}m" BaseRecalibrator \
            --reference "~{reference.fasta}" \
            --input "~{inputBam}" \
            --use-original-qualities \
            --output ~{outputRecalibrationReport} \
            --known-sites ~{dbsnp.file} \
            $(printf ' --known-sites %s ' $(printf '%s\n' ${KNOWNSITES[@]}))
            
            #-L Interval list separated format
    >>>
    
    output {
        File recalibrationReport = outputRecalibrationReport
    }

    runtime {
        memory: select_first([memoryGb * 1024,6*1024])
        timeMinutes: timeMinutes
    }
}

#gatherbsqr reports

task ApplyBQSR {
  input {
    File inputBam
    File inputBai
    String gatkModule = "GATK"
    String outputBamBasename
    File recalibrationReport
    Reference reference
    Int? memoryGb = "4"
    Int? javaXmxMemoryMb = floor(memoryGb*0.9*1024)
    Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
    Int disk = ceil(size(inputBam, "M")*2.6)
  }
  #Float referenceSize = size(reference.fasta, "GiB") + size(reference.dict, "GiB") + size(reference.fai, "GiB")
  #Int DiskSize = ceil((size(inputBam, "GiB") * 3 ) + referenceSize) + 20
  command {
    ml ~{gatkModule}
    gatk --java-options "-XX:+PrintFlagsFinal \
        -XX:+PrintGCDetails -Xloggc:gc_log.log \
        -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m -Xmx~{javaXmxMemoryMb}m" \
        ApplyBQSR \
        --create-output-bam-md5 true\
        --create-output-bam-index true\
        --add-output-sam-program-record \
        -R ~{reference.fasta} \
        -I ~{inputBam} \
        --use-original-qualities \
        -O ~{outputBamBasename}.bam \
        -bqsr ~{recalibrationReport} 
  }
  runtime {
    memory: memoryGb
    timeMinutes: timeMinutes
  }
  output {
    File bam = outputBamBasename + ".bam"
    File bai = outputBamBasename + ".bai"
    File bammd5sum = outputBamBasename + ".bam.md5"
  }
}

task HaplotypeCallerGVcf {
    input {
        IndexedFile inputBam
        #File inputBai
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
        Int timeMinutes = 1 + ceil(size(inputBam.file, "G")) * 120
        Int? javaXmxMemoryMb = floor(memoryGb*0.9*1024)
        Int disk = ceil(size(inputBam.file, "M")*1.2)
    }

    String vcfSuffix =  if makeGvcf then ".g.vcf" else ".vcf"
    String bamoutArg =  if makeBamOut then "-bamout " + outputVcfBasename + ".bamout.bam" else ""
    File bamIn = inputBam.file
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
        set -e

        ml ~{gatkModule}
        gatk --java-options "-Xmx${javaXmxMemoryMb}m -Xms${javaXmxMemoryMb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            HaplotypeCaller \
            -R ~{reference.fasta} \
            -I ~{bamIn} \
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
        File vcf = outputVcfBasename + vcfSuffix
        File vcfIdx = outputVcfBasename + vcfSuffix + ".idx"
        IndexedFile vcfOut = {
          "file" : outputVcfBasename + vcfSuffix,
          "index" : outputVcfBasename + vcfSuffix + ".idx"
        }
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
    }
}

#add genomicsdb import#

task GenotypeGVcf {
    input{
        Array[File] inputGVcfs
        Reference reference
        String gatkModule = "GATK"
        Int memoryGb = "4"
        Int javaMemoryGb = memoryGb - 1
        Int timeMinutes = 1 + ceil(size(inputGVcfs, "G")) * 120

      }
    command <<<
        ml ~{gatkModule}
      

    >>>
    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
    }

}
#SPLIT_TO_N_READS
#https://github.com/broadinstitute/warp/blob/develop/tasks/broad/Alignment.wdl#L128