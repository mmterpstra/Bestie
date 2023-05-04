version 1.0

task FastqToUnmappedBam {
    input {
        File inputFastq1
        File? inputFastq2
        Int? memoryGb = "16"
        String picardModule = "picard"
        String sampleName = "test"
        String readgroup = "A"
        String platform = "ILLUMINA"
        String platformUnit = "run_barcode.lane"
        String readGroupName = "flowcell_run_barcode.lane"
        #String? platformModel = "NextSeq?"    
        String outputUnalignedBam = "unaligned_test.sam"
        Int timeMinutes = 1 + ceil(size(inputFastq1, "G")) * 40
    }
    command {
        set -e
        module load ~{picardModule} && \
        java -Xmx1g \
            -jar $EBROOTPICARD/picard.jar FastqToSam \
            FASTQ=~{inputFastq1} \
            ~{"FASTQ2=" + inputFastq2} \
            SAMPLE_NAME=~{sampleName} \
            PLATFORM=~{platform} \
            RUN_DATE="$(date --rfc-3339=date)" \
            PLATFORM_UNIT=~{platformUnit} \
            READ_GROUP_NAME=~{readGroupName} \
            OUTPUT=~{outputUnalignedBam}
    }

    output {
        File unalignedBam = outputUnalignedBam
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
    }
}

task SortSam {
    input {
        File inputBam
        Int? memoryGb = "5"
        String picardModule = "picard"
        String outputBamBasename
        #String? platformModel = "NextSeq?"    
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 40
    }
    command {
        set -e
        module load ~{picardModule} && \
        java -Xmx4g \
            -jar $EBROOTPICARD/picard.jar SortSam \
            INPUT=~{inputBam} \
            OUTPUT=~{outputBamBasename}.bam \
            SORT_ORDER="coordinate" \
            CREATE_INDEX=true \
            CREATE_MD5_FILE=true \
            MAX_RECORDS_IN_RAM=300000
    }

    output {
        File bam = outputBamBasename + ".bam"
        File bai = outputBamBasename + ".bai"
        File md5 = outputBamBasename + ".bam.md5"
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
    }
}

task MarkDuplicates {
    input {
        Array[File] inputBams
        String outputBamBasename
        String outputMetrics
        String picardModule = "picard"
        Int? memoryGb = "5"
        Int timeMinutes = 1 + ceil(size(inputBams, "G")) * 40
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
      ml ~{picardModule}
      java -Xmx4g -jar $EBROOTPICARD/picard.jar \
      MarkDuplicates \
      INPUT=~{sep=' INPUT=' inputBams} \
      OUTPUT=~{outputBamBasename}.bam \
      METRICS_FILE=~{outputMetrics} \
      VALIDATION_STRINGENCY=SILENT \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER="queryname" \
      CLEAR_DT="false" \
      ADD_PG_TAG_TO_READS=false
    }
    
    output {
        File bam = outputBamBasename + ".bam"
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
    }
}

#SPLIT_TO_N_READS
#https://github.com/broadinstitute/warp/blob/develop/tasks/broad/Alignment.wdl#L128