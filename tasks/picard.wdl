version 1.0

import "../structs.wdl" as structs


task FastqToUnmappedBam {
    input {
        File inputFastq1
        File? inputFastq2
        Int memoryGb = "1"
        Int compressionLevel = 5
        Int javaXmxMemoryMb = floor((memoryGb-0.5)*0.95*1024)
        String picardModule = "picard"
        String sampleName = "test"
        String libraryName
        String readgroup = "A"
        String platform = "ILLUMINA"
        String platformUnit = "run_barcode.lane"
        String readGroupName = "flowcell_run_barcode.lane"
        #String? platformModel = "NextSeq?"    
        String outputUnalignedBam = "unaligned_test.sam"
        Int timeMinutes = 1 + ceil(size(inputFastq1, "G")) * 40
        #Possible values: {unsorted, queryname, coordinate, duplicate, unknown} #
        String sortOrder = "queryname" 
    }
    command {
        set -e
        module load ~{picardModule} && \
        java -Xmx~{javaXmxMemoryMb}m \
            -jar $EBROOTPICARD/picard.jar FastqToSam \
            FASTQ=~{inputFastq1} \
            ~{"FASTQ2=" + inputFastq2} \
            SAMPLE_NAME=~{sampleName} \
            ~{"LIBRARY_NAME=" + libraryName}\
            PLATFORM=~{platform} \
            RUN_DATE="$(date --rfc-3339=date)" \
            PLATFORM_UNIT=~{platformUnit} \
            READ_GROUP_NAME=~{readGroupName} \
            OUTPUT=~{outputUnalignedBam} \
            SORT_ORDER=~{sortOrder} \
            COMPRESSION_LEVEL=~{compressionLevel}

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
        Int memoryGb = "2"
        Int compressionLevel = 5
        Int javaXmxMemoryMb = floor((memoryGb-0.5)*0.95*1024)
        String picardModule = "picard"
        Int disk = ceil(size(inputBam, "M")*2.1)
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 40
    }
    command {
        set -e
        module load ~{picardModule} 
        mkdir -p ~{outputFastqDirBase}
        java -Xmx~{javaXmxMemoryMb}m \
            -jar $EBROOTPICARD/picard.jar SamToFastq \
            INPUT=~{inputBam} \
            OUTPUT_PER_RG=true \
            COMPRESS_OUTPUTS_PER_RG=true \
            OUTPUT_DIR=~{outputFastqDirBase} \
            COMPRESSION_LEVEL=~{compressionLevel}

    }
    output {
        File fastq1gz = select_first(glob(outputFastqDirBase + "/*_1.fastq.gz"))
        File? fastq2gz = select_first(glob(outputFastqDirBase + "/*_2.fastq.gz"))
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
        Int memoryGb = "15"
        Int compressionLevel = 5
        Int javaXmxMemoryMb = floor((memoryGb-0.5)*0.95*1024)
        String picardModule = "picard"
        String outputBamBasename
        #String? platformModel = "NextSeq?"    
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 40
        Int disk = ceil(size(inputBam, "M")*2.1)
        #Possible values: {unsorted, queryname, coordinate, duplicate, unknown} #
        String sortOrder = "coordinate" 
    }
    Boolean createIndex = if sortOrder=="coordinate" then true else false
    command {
        set -e
        module load ~{picardModule} && \
        java -Xmx~{javaXmxMemoryMb}m \
            -jar $EBROOTPICARD/picard.jar SortSam \
            INPUT=~{inputBam} \
            OUTPUT=~{outputBamBasename}.bam \
            SORT_ORDER=~{sortOrder} \
            ~{true="CREATE_INDEX=true " false="" createIndex} \
            CREATE_MD5_FILE=true \
            MAX_RECORDS_IN_RAM=300000 \
            COMPRESSION_LEVEL=~{compressionLevel}

    }

    output {
        File bam = outputBamBasename + ".bam"
        File? bai = outputBamBasename + ".bai"
        File? md5 = outputBamBasename + ".bam.md5"
    }

    runtime {
        memory: select_first([memoryGb * 1024, 4*1024])
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
        Int memoryGb = 16
        Int compressionLevel = 5
        Int javaXmxMemoryMb = floor((memoryGb-0.5)*0.95*1024)
        Int timeMinutes = 1 + ceil(size(inputBams, "G")) * 60
        Int disk = ceil(size(inputBams, "M")*1.2)
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
      ml ~{picardModule}
      java -Xmx~{javaXmxMemoryMb}m -jar $EBROOTPICARD/picard.jar \
      MarkDuplicates \
      INPUT=~{sep=' INPUT=' inputBams} \
      OUTPUT=~{outputBamBasename}.bam \
      METRICS_FILE=~{outputMetrics} \
      VALIDATION_STRINGENCY=SILENT \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER="queryname" \
      CLEAR_DT="false" \
      ADD_PG_TAG_TO_READS=false \
      COMPRESSION_LEVEL=~{compressionLevel}
    }
    
    output {
        File bam = outputBamBasename + ".bam"
        File metrics = outputMetrics
    }

    runtime {
        memory: select_first([memoryGb * 1024, 4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}
task UmiAwareMarkDuplicatesWithMateCigar {
    input {
        Array[File] inputBams
        String outputBamBasename
        String outputMetrics
        String outputUMIMetrics
        String UmiTagName="RX"
        String picardModule = "picard"
        Int memoryGb = 16
        Int compressionLevel = 5
        Int javaXmxMemoryMb = floor((memoryGb-0.5)*0.95*1024)
        Int timeMinutes = 1 + ceil(size(inputBams, "G")) * 120
        Int disk = ceil(size(inputBams, "M")*1.2)
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {

        set -e -o pipefail
        
        ml ~{picardModule}

        mkfifo ~{outputBamBasename}_umimarkdupmatecigar.fifo.bam

        java -Xmx1078m -jar $EBROOTPICARD/picard.jar MergeSamFiles \
            INPUT=~{sep=' INPUT=' inputBams} \
            OUTPUT=/dev/stdout \
            COMPRESSION_LEVEL=0 | \
        java -Xmx4078m -jar $EBROOTPICARD/picard.jar SortSam \
            INPUT=/dev/stdin \
            OUTPUT=~{outputBamBasename}_umimarkdupmatecigar.fifo.bam \
            SORT_ORDER=coordinate \
            MAX_RECORDS_IN_RAM=300000 \
            COMPRESSION_LEVEL=0 &
        
        java -Xmx~{javaXmxMemoryMb}m -jar $EBROOTPICARD/picard.jar \
            UmiAwareMarkDuplicatesWithMateCigar \
            INPUT=~{outputBamBasename}_umimarkdupmatecigar.fifo.bam \
            OUTPUT=~{outputBamBasename}.bam \
            METRICS_FILE=~{outputMetrics} \
            UMI_TAG_NAME=~{UmiTagName} \
            UMI_METRICS_FILE=~{outputUMIMetrics}\
            VALIDATION_STRINGENCY=SILENT \
            CLEAR_DT="false" \
            CREATE_INDEX=true \
            COMPRESSION_LEVEL=~{compressionLevel}

        rm -v ~{outputBamBasename}_umimarkdupmatecigar.fifo.bam
    }
    
    output {
        File bam = outputBamBasename + ".bam"
        File bai = outputBamBasename + ".bai"
        File metrics = outputMetrics
        File umiMetrics = outputUMIMetrics
    }

    runtime {
        memory: select_first([memoryGb * 1024 + 6000, 4*1024 + 6000])
        timeMinutes: timeMinutes
        disk: disk
    }
}
task MergeSamFiles {
    input {
        Array[File] inputBams
        String outputBamBasename
        String picardModule = "picard"
        Int memoryGb = "5"
        Int compressionLevel = 5
        Int javaXmxMemoryMb = floor((memoryGb-0.5)*0.95*1024)
        Int timeMinutes = 1 + ceil(size(inputBams, "G")) * 40
        Int disk = 1 + ceil(size(inputBams, "G") * 2.1)
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
        set -e -o pipefail

        ml ~{picardModule}
        java -Xmx~{javaXmxMemoryMb}m -XX:ParallelGCThreads=4 -jar $EBROOTPICARD/picard.jar MergeSamFiles \
            INPUT=~{sep=' INPUT=' inputBams} \
            SORT_ORDER=coordinate \
            CREATE_INDEX=true \
            USE_THREADING=true \
            TMP_DIR=./ \
            MAX_RECORDS_IN_RAM=6000000 \
            COMPRESSION_LEVEL=~{compressionLevel} \
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
        Int memoryGb = "5"
        Int javaXmxMemoryMb = floor((memoryGb-0.5)*0.95*1024)
        Int timeMinutes = 1 + ceil(size(inputIntervalListFile, "G")) * 20
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
        ml ~{picardModule}
        #mkdir -p ~{outputPrefix}

        java -Xmx~{javaXmxMemoryMb}m -jar $EBROOTPICARD/picard.jar \
        IntervalListTools \
        INPUT="~{inputIntervalListFile}" \
        OUTPUT="padded_list.interval_list" \
        PADDING=~{padding} \
        SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
        UNIQUE=true \
        COMMENT="Added padding of ~{padding} bp and merge overlapping and adjacent intervals to create a list of unique intervals PADDING=~{padding} UNIQUE=true"

        mkdir -p ~{outputPrefix}_scatter
        mkdir -p scatter_list
        java -Xmx~{javaXmxMemoryMb}m -jar $EBROOTPICARD/picard.jar \
        IntervalListTools \
        INPUT="~{inputIntervalListFile}" \
        OUTPUT="scatter_list" \
        PADDING=~{padding} \
        SCATTER_COUNT=~{targetScatter} \
        SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
        UNIQUE=true \
        COMMENT="Added padding of ~{padding} bp and merge overlapping and adjacent intervals to create a list of unique intervals PADDING=~{padding} UNIQUE=true"  
        FILEINDEX=0
        for FILE in ./scatter_list/*/*.interval_list; do
            mv $FILE "$(dirname $FILE)""/""$FILEINDEX""_""$(basename $FILE)"
            FILEINDEX=$((FILEINDEX+1))
        done
    }
    
    output {
        File paddedIntervalList = "padded_list.interval_list"
        #Array[File] paddedScatteredIntervalList = glob(outputPrefix + "_scatter/temp_*_of_"+ targetScatter +"/scattered.interval_list")
        Array[File] paddedScatteredIntervalList = glob("scatter_list/*/*.interval_list")
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
    }
}

task GatherVcfs {
    input {
        Array[File] inputVcfs
        Boolean createIndex = true
        String outputPrefix
        String vcfSuffix = ".vcf.gz"
        String picardModule = "picard"
        Int memoryGb = "5"
        Int javaXmxMemoryMb = floor((memoryGb-0.5)*0.95*1024)
        Int timeMinutes = 1 + ceil(size(inputVcfs, "G")) * 60
    }
    command {
        module load ~{picardModule}
        java -Xmx~{javaXmxMemoryMb}m -jar $EBROOTPICARD/picard.jar \
        GatherVcfs \
        ~{true="CREATE_INDEX=true " false="" createIndex} \
        INPUT=~{sep=' INPUT=' inputVcfs} \
        OUTPUT=~{outputPrefix}~{vcfSuffix} \
    }
    output {
        File outputVcf = outputPrefix + vcfSuffix
        #File outputVcfIdx = outputPrefix + vcfSuffix+ ".tbi"
        #IndexedFile vcfOut = {
        #  "file" : outputPrefix + vcfSuffix,
        #  "index" : outputPrefix + vcfSuffix + ".tbi"
        #}
    }
    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
    }
}

task GatherVcfsIndexed {
    input {
        Array[File] inputVcfs
        Boolean createIndex = true
        String outputPrefix
        String vcfSuffix = ".vcf.gz"
        String picardModule = "picard"
        Int memoryGb = "5"
        Int javaXmxMemoryMb = floor((memoryGb-0.5)*0.95*1024)
        Int timeMinutes = 1 + ceil(size(inputVcfs, "G")) * 60
    }
    command {
        set -o pipefail
        module load ~{picardModule}
        java -Xmx~{javaXmxMemoryMb}m -jar $EBROOTPICARD/picard.jar \
        GatherVcfs \
        INPUT=$(realpath ~{sep=') INPUT=$(realpath ' inputVcfs}) \
        OUTPUT=/dev/stdout | \
        java -Xmx~{javaXmxMemoryMb}m -jar $EBROOTPICARD/picard.jar \
        SortVcf \
        INPUT=/dev/stdin \
        ~{true="CREATE_INDEX=true " false="" createIndex} \
        OUTPUT=~{outputPrefix}~{vcfSuffix}

    }
    output {
        File outputVcf = outputPrefix + vcfSuffix
        File outputVcfIdx = outputPrefix + vcfSuffix+ ".tbi"
        IndexedFile vcfOut = {
          "file" : outputPrefix + vcfSuffix,
          "index" : outputPrefix + vcfSuffix + ".tbi"
        }
    }
    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
    }
}
task SortVcfsIndexed {
    input {
        Array[File] inputVcfs
        Boolean createIndex = true
        String outputPrefix
        String vcfSuffix = ".vcf.gz"
        String picardModule = "picard"
        Int memoryGb = "5"
        Int javaXmxMemoryMb = floor((memoryGb-0.5)*0.95*1024)
        Int timeMinutes = 1 + ceil(size(inputVcfs, "G")) * 60
    }
    command {
        set -o pipefail
        module load ~{picardModule}
        java -Xmx~{javaXmxMemoryMb}m -jar $EBROOTPICARD/picard.jar \
        SortVcf \
        INPUT= ~{sep=' INPUT= ' inputVcfs} \
        ~{true="CREATE_INDEX=true " false="" createIndex} \
        OUTPUT=~{outputPrefix}~{vcfSuffix}

    }
    output {
        File outputVcf = outputPrefix + vcfSuffix
        File outputVcfIdx = outputPrefix + vcfSuffix+ ".tbi"
        IndexedFile vcfOut = {
          "file" : outputPrefix + vcfSuffix,
          "index" : outputPrefix + vcfSuffix + ".tbi"
        }
    }
    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
    }
}


task SplitSamByNumberOfReads {
    input {
        File inputBam
        Int memoryGb = "4"
        Int compressionLevel = 5
        Int javaXmxMemoryMb = floor((memoryGb-0.5)*0.95*1024)
        String picardModule = "picard"
        String samtoolsModule = "SAMtools"
        String outputBamBaseDir
        #50m or bigger
        Int numberOfReads = 5000000
        #String? platformModel = "NextSeq?"    
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
        Int disk = ceil(size(inputBam, "M")*2.1)
        #Possible values: {unsorted, queryname, coordinate, duplicate, unknown} #
        String sortOrder = "coordinate" 
    }
    command {
        set -e
        TOTALREADS=$(module load ~{samtoolsModule} && samtools view -c  ~{inputBam})

        module load ~{picardModule}&& \
        mkdir -p ~{outputBamBaseDir} && \
        java -Xmx~{javaXmxMemoryMb}m \
            -jar $EBROOTPICARD/picard.jar SplitSamByNumberOfReads \
            INPUT=~{inputBam} \
            OUTPUT=~{outputBamBaseDir} \
            SPLIT_TO_N_READS=~{numberOfReads} \
            COMPRESSION_LEVEL=~{compressionLevel} \
            TOTAL_READS_IN_INPUT=$TOTALREADS
        
    }

    output {
        Array[File] bams = glob(outputBamBaseDir+"/*.bam")
    }

    runtime {
        memory: select_first([memoryGb * 1024, 4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}