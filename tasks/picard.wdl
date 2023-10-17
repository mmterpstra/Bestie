version 1.0

task FastqToUnmappedBam {
    input {
        File inputFastq1
        File? inputFastq2
        Int? memoryGb = "2"
        String picardModule = "picard"
        String sampleName = "test"
        String readgroup = "A"
        String platform = "ILLUMINA"
        String platformUnit = "run_barcode.lane"
        String readGroupName = "flowcell_run_barcode.lane"
        #String? platformModel = "NextSeq?"    
        String outputUnalignedBam = "unaligned_test.sam"
        Int timeMinutes = 1 + ceil(size(inputFastq1, "G")) * 120
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

task SamToFastq {
    input{
        File inputBam
        String outputFastqDirBase
        Int? memoryGb = "2"
        String picardModule = "picard"
        Int disk = ceil(size(inputBam, "M")*2.1)
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 60
    }
    command {
        set -e
        module load ~{picardModule} && \
        java -Xmx1g \
            -jar $EBROOTPICARD/picard.jar FastqToSam \
            INPUT=~{inputBam} \
            OUTPUT_PER_RG=true \
            OUTPUTDIR=~{outputFastqDirBase}
    }
    output {
        File fastq1gz = select_first(glob(outputFastqDirBase + "/*R1.fastq.gz"))
        File? fastq2gz = select_first(glob(outputFastqDirBase + "/*R2.fastq.gz"))
    }
    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}

task SortSam {
    input {
        File inputBam
        Int? memoryGb = "5"
        String picardModule = "picard"
        String outputBamBasename
        #String? platformModel = "NextSeq?"    
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
        Int disk = ceil(size(inputBam, "M")*2.1)
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
        disk: disk
    }
}

task MarkDuplicates {
    input {
        Array[File] inputBams
        String outputBamBasename
        String outputMetrics
        String picardModule = "picard"
        Int? memoryGb = "5"
        Int timeMinutes = 1 + ceil(size(inputBams, "G")) * 120
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
        File metrics = outputMetrics
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
    }
}

task MergeSamFiles {
    input {
        Array[File] inputBams
        String outputBamBasename
        String picardModule = "picard"
        Int? memoryGb = "5"
        Int timeMinutes = 1 + ceil(size(inputBams, "G")) * 120
        Int disk = 1 + ceil(size(inputBams, "G") * 2.1)
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
        ml ~{picardModule}
        java -Xmx4g -XX:ParallelGCThreads=4 -jar $EBROOTPICARD/picard.jar MergeSamFiles \
            INPUT=~{sep=' INPUT=' inputBams} \
            SORT_ORDER=coordinate \
            CREATE_INDEX=true \
            USE_THREADING=true \
            TMP_DIR=./ \
            MAX_RECORDS_IN_RAM=6000000 \
            OUTPUT=${outputBamBasename}.bam
    }
    
    output {
        File bam = outputBamBasename + ".bam"
        File bai = outputBamBasename + ".bai"
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}

task SplitAndPadIntervals {
    input {
        File inputIntervalListFile
        String outputPrefix
        Int padding = 150
        Int targetScatter = 50
        String picardModule = "picard"
        Int? memoryGb = "4"
        Int timeMinutes = 1 + ceil(size(inputIntervalListFile, "G")) * 120
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
        ml ~{picardModule}
        #mkdir -p ~{outputPrefix}

        java -Xmx4g -jar $EBROOTPICARD/picard.jar \
        IntervalListTools \
        INPUT="~{inputIntervalListFile}" \
        OUTPUT="~{outputPrefix}.interval_list" \
        PADDING=~{padding} \
        SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
        UNIQUE=true \
        COMMENT="Added padding of ~{padding} bp and merge overlapping and adjacent intervals to create a list of unique intervals PADDING=~{padding} UNIQUE=true"

        mkdir -p ~{outputPrefix}_scatter
        
        java -Xmx4g -jar $EBROOTPICARD/picard.jar \
        IntervalListTools \
        INPUT="~{inputIntervalListFile}" \
        OUTPUT="~{outputPrefix}_scatter" \
        PADDING=~{padding} \
        SCATTER_COUNT=~{targetScatter} \
        SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
        UNIQUE=true \
        COMMENT="Added padding of ~{padding} bp and merge overlapping and adjacent intervals to create a list of unique intervals PADDING=~{padding} UNIQUE=true"  
    }
    
    output {
        File paddedIntervalList = outputPrefix + ".interval_list"
        #Array[File] =  outputPrefix + ".interval_list"
        Array[File] paddedScatteredIntervalList = glob(outputPrefix + "_scatter/temp_*_of_"+ targetScatter +"/scattered.interval_list")
        #Int interval_count = read_int(stdout())
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
    }
}
