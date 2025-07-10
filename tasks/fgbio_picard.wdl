version 1.0

import "../structs.wdl" as structs

task FastqToUnmappedBamPicardSorted {
    input {
        File inputFastq1
        File? inputFastq2
        File? inputUmiFastq1
        Int memoryGb = "4"
        Int javaPicardXmxMemoryMb = floor(2*0.95*1024)
        Int javaFgbioXmxMemoryMb = floor(1*0.95*1024)
        #needs fgbio2
        String fgbioModule = "fgbio/2.2.1-Java-8-LTS"
        #picard2
        String picardModule = "picard"
        String sampleName = "test"
        String readgroup = "A"
        String library = "sample_barcode1+barcode2"
        String platform = "ILLUMINA"
        String platformUnit = "run_barcode.lane"
        String readGroupName = "flowcell_run_barcode.lane"
        #String? platformModel = "NextSeq?"    
        String outputUnalignedBam = "unaligned_test.sam"
        Int timeMinutes = 1 + ceil(size(inputFastq1, "G")) * 120
        #default +T or "5M2S+T" for twist datasets
        String? readStructureFastq1 = "+T"
        #default +T or "5M2S+T" for twist datasets
        String? readStructureFastq2 = "+T"
        #default +M when defined
        String? readStructureFastqUmi = ""
        Boolean? extractUmisFromReadNames = false
        Boolean? sort = false    
    }
    command {
        set -e -o pipefail
        
        
        READSTRUCTURES="~{readStructureFastq1} ~{readStructureFastq2} ~{readStructureFastqUmi}"
        UMIQUALTAG=""
        if [[ $READSTRUCTURES == *M* ]]; then
            if ~{true="true" false="false" extractUmisFromReadNames} ; then
                1>&2 echo "## ERROR ## unspecified behavior: 'extract umitags from readnames' clashing with 'umi spec in readstructures'. "
                exit 1
            fi
            1>&2 echo  "The readstructures contain one or more 'M'. tags. Setting umi qual tag."
            UMIQUALTAG="--umi-qual-tag RQ "
        fi
        

        (
            module load ~{fgbioModule} && java -Xmx~{javaFgbioXmxMemoryMb}m \
            --compression 0 \
            -jar $EBROOTFGBIO/lib/fgbio-$(echo ~{fgbioModule} | perl -wpe 's/fgbio\/([\d.]+).*/$1/g').jar FastqToBam \
            --input ~{inputFastq1} ~{inputFastq2} ~{inputUmiFastq1} \
            --read-structures ~{readStructureFastq1} ~{readStructureFastq2} ~{readStructureFastqUmi} \
            ~{true=" --extract-umis-from-read-names " false="" extractUmisFromReadNames} \
            --umi-tag RX \
            $UMIQUALTAG \
            --sample ~{sampleName} \
            --library ~{library} \
            --platform ~{platform} \
            --run-date "$(date --rfc-3339=date)" \
            --platform-unit ~{platformUnit} \
            --output /dev/stdout 
        ) | (
            module load ~{picardModule} && java -Xmx~{javaPicardXmxMemoryMb}m \
                -jar $EBROOTPICARD/picard.jar SortSam \
                INPUT=/dev/stdin \
                OUTPUT=~{outputUnalignedBam} \
                SORT_ORDER=queryname \
                MAX_RECORDS_IN_RAM=300000
        )
    }

    output {
        File ubam = outputUnalignedBam
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
    }
}

task FastqToFgUnmappedBamPicardSortedScattered {
    input {
        File inputFastq1
        File? inputFastq2
        File? inputUmiFastq1
        Int memoryGb = "4"
        Int javaPicardSortXmxMemoryMb = floor(2*0.95*1024)
        Int javaPicardSplitXmxMemoryMb = floor(2*0.95*1024)
        Int javaFgbioXmxMemoryMb = floor(2*0.95*1024)
        #needs fgbio2
        String fgbioModule = "fgbio/2.2.1-Java-8-LTS"
        #picard2
        String picardModule = "picard"
        String sampleName = "test"
        String readgroup = "A"
        String library = "sample_barcode1+barcode2"
        String platform = "ILLUMINA"
        String platformUnit = "run_barcode.lane"
        String readGroupName = "flowcell_run_barcode.lane"
        #String? platformModel = "NextSeq?"    
        #String outputUnalignedBam = "unaligned_test.sam"
        Int timeMinutes = 1 + ceil(size(inputFastq1, "G")) * 120
        #default +T or "5M2S+T" for twist datasets
        String? readStructureFastq1 = "+T"
        #default +T or "5M2S+T" for twist datasets
        String? readStructureFastq2 = "+T"
        #default +M when defined
        String? readStructureFastqUmi = ""
        Boolean? extractUmisFromReadNames = false
        Int numberOfReads = 50000000
        String outputBamBaseDir
    }
    command {
        set -e -o pipefail

        TOTALREADS=$(gzip -qdc ~{inputFastq1} ~{inputFastq2} | \
            wc -l | \
            python3 -c 'import sys; print("\n".join(str(int(x)//4) for x in sys.stdin))' )
        SPLIT_TO_N_READS=~{numberOfReads}
        if [ "~{numberOfReads}" -gt "$TOTALREADS" ]; then

            SPLIT_TO_N_READS=$TOTALREADS
        fi

        READSTRUCTURES="~{readStructureFastq1} ~{readStructureFastq2} ~{readStructureFastqUmi}"
        UMIQUALTAG=""
        if [[ $READSTRUCTURES == *M* ]]; then
            if ~{true="true" false="false" extractUmisFromReadNames} ; then
                1>&2 echo "## ERROR ## unspecified behavior: 'extract umitags from readnames' clashing with 'umi spec in readstructures'. "
                exit 1
            fi
            1>&2 echo  "The readstructures contain one or more 'M'. tags. Setting umi qual tag."
            UMIQUALTAG="--umi-qual-tag RQ "
        fi
        mkdir -p "~{outputBamBaseDir}"

        (
            module load ~{fgbioModule} && java -Xmx~{javaFgbioXmxMemoryMb}m \
            -jar $EBROOTFGBIO/lib/fgbio-$(echo ~{fgbioModule} | perl -wpe 's/fgbio\/([\d.]+).*/$1/g').jar \
            --compression 0 \
            FastqToBam \
            --input ~{inputFastq1} ~{inputFastq2} ~{inputUmiFastq1} \
            --read-structures ~{readStructureFastq1} ~{readStructureFastq2} ~{readStructureFastqUmi} \
            ~{true=" --extract-umis-from-read-names " false="" extractUmisFromReadNames} \
            --umi-tag RX \
            $UMIQUALTAG \
            --sample ~{sampleName} \
            --library ~{library} \
            --platform ~{platform} \
            --run-date "$(date --rfc-3339=date)" \
            --platform-unit ~{platformUnit} \
            --output /dev/stdout 
        ) | (
            module load ~{picardModule} && java -Xmx~{javaPicardSortXmxMemoryMb}m \
                -jar $EBROOTPICARD/picard.jar SortSam \
                INPUT=/dev/stdin \
                OUTPUT=/dev/stdout \
                COMPRESSION_LEVEL=0 \
                SORT_ORDER=queryname \
                MAX_RECORDS_IN_RAM=300000
        ) | ( 
            module load ~{picardModule} && \
            java -Xmx~{javaPicardSplitXmxMemoryMb}m \
                -jar $EBROOTPICARD/picard.jar SplitSamByNumberOfReads \
                INPUT=/dev/stdin \
                OUTPUT="~{outputBamBaseDir}" \
                SPLIT_TO_N_READS=$SPLIT_TO_N_READS \
                TOTAL_READS_IN_INPUT=$TOTALREADS
        )
    }

    output {
        Array[File] ubams = glob(outputBamBaseDir+"/*.bam")
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
    }
}
