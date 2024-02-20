version 1.0

task FastqToUnmappedBam {
    input {
        File inputFastq1
        File? inputFastq2
        File? inputUmiFastq1
        String readStr
        Int memoryGb = "1"
        Int javaXmxMemoryMb = floor(memoryGb*0.95*1024)
        String fgbioModule = "fgbio"
        String sampleName = "test"
        String readgroup = "A"
        String library = "sample_barcode1+barcode2"
        String platform = "ILLUMINA"
        String platformUnit = "run_barcode.lane"
        String readGroupName = "flowcell_run_barcode.lane"
        #String? platformModel = "NextSeq?"    
        String outputUnalignedBam = "unaligned_test.sam"
        Int timeMinutes = 1 + ceil(size(inputFastq1, "G")) * 120
    }
    command {
        set -e
        module load ~{fgbioModule} && \
        java -Xmx~{javaXmxMemoryMb}m \
            -jar $EBROOTFGBIO/lib/fgbio-$(echo ~{fgbioModule} | perl -wpe 's/fgbio\/([\d.]+).*/$1/g').jar FastqToBam \
            --input ~{inputFastq1} ~{inputFastq2} ~{inputUmiFastq1} \
            --read-structures +T +T +M \
            --umi-tag RX \
            --umi-qual-tag RQ \
            --sample ~{sampleName} \
            --library ~{library} \
            --platform ~{platform} \
            --run-date "$(date --rfc-3339=date)" \
            --platform-unit ~{platformUnit} \
            --output ~{outputUnalignedBam}
    }

    output {
        File unalignedBam = outputUnalignedBam
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
    }
}

task ExtractUmisFromBam {
    input {
        File inputBam
        String outputBamBasename
        String fgbioModule = "fgbio"
        Int memoryGb = "4"
        Int javaMemoryMb = ceil(memoryGb * 1024 * 0.85)
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
        Int disk = ceil(size(inputBam, "M")*2.1)
    }
    command {
        set -e
        set -o pipefail
        module load ~{fgbioModule} && \
        java -Xmx~{javaMemoryMb}m \
            -jar $EBROOTFGBIO/lib/fgbio-$(echo ~{fgbioModule} | perl -wpe 's/fgbio\/([\d.]+).*/$1/g').jar \
            ExtractUmisFromBam \
            --input=~{inputBam} \
            --output=~{outputBamBasename}.bam \
            --read-structure="5M2S+T" \
            --read-structure="5M2S+T" \
            --molecular-index-tags=ZA \
            --molecular-index-tags=ZB \
            --single-tag=RX
    }

    output {
        File bam = outputBamBasename + ".bam"
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}
task GroupReadsByUmi {
    input {
        File inputBam
        String outputBamBasename
        String fgbioModule = "fgbio"
        String groupUmiReadsStrategy = "paired" #can also be identity, edit or adjacency. See https://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html
        Int memoryGb = "4"
        Int javaMemoryMb = ceil(memoryGb * 1024 * 0.85)
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
        Int disk = ceil(size(inputBam, "M")*2.1)
    }
    command {
        set -e
        module load ~{fgbioModule} && \
        java -Xmx~{javaMemoryMb}m -jar $EBROOTFGBIO/lib/fgbio-$(echo ~{fgbioModule} | perl -wpe 's/fgbio\/([\d.]+).*/$1/g').jar GroupReadsByUmi \
            --strategy=~{groupUmiReadsStrategy} \
            --input=~{inputBam} \
            --output=~{outputBamBasename}".bam" \
            --family-size-histogram=~{outputBamBasename}".hist" \
            --raw-tag=RX \
            --min-map-q=10 \
            --edits=1
    }

    output {
        File bam = outputBamBasename + ".bam"
        File familySizeHistogram = outputBamBasename + ".hist"
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}

task CallDuplexConsensusReads {
    input {
        File inputBam
        String outputBamBasename
        String fgbioModule = "fgbio"
        Int memoryGb = "5"
        Int javaMemoryMb = ceil(memoryGb * 1024 * 0.85)
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
        Int disk = ceil(size(inputBam, "M")*2.1)
    }
    command {
        set -e
        module load ~{fgbioModule} && \
        java -Xmx~{javaMemoryMb}m -jar $EBROOTFGBIO/lib/fgbio-$(echo ~{fgbioModule} | perl -wpe 's/fgbio\/([\d.]+).*/$1/g').jar CallDuplexConsensusReads \
            --input=~{inputBam} \
            --output="~{outputBamBasename}.bam" \
            --error-rate-pre-umi=45 \
            --error-rate-post-umi=30 \
            --min-input-base-quality=30 \
            --min-reads 2 1 1
    }

    output {
        File bam = outputBamBasename + ".bam"
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}
#wip
task FilterDuplexConsensusReads {
    input {
        File inputBam
        String outputBamBasename
        String fgbioModule = "fgbio"
        Int memoryGb = "5"
        Int javaMemoryMb = ceil(memoryGb * 1024 * 0.85)
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
        Int disk = ceil(size(inputBam, "M")*2.1)
    }
    command {
        set -e
        module load ~{fgbioModule} && \
        java -Xmx~{javaMemoryMb}m -jar $EBROOTFGBIO/lib/fgbio-$(echo ~{fgbioModule} | perl -wpe 's/fgbio\/([\d.]+).*/$1/g').jar FilterDuplexConsensusReads \
            --input=~{inputBam} \
            --output="~{outputBamBasename}.bam" 
            #\
            #    --error-rate-pre-umi=45 \
            #    --error-rate-post-umi=30 \
            #    --min-input-base-quality=30 \
            #    --min-reads 2 1 1
            #--reverse-per-base-tags=true \
            #--min-reads=3 \
            #-max-read-error-rate 0.05 \
            #-min-base-quality 40 \
            #-max-base-error-rate 0.1 \
            #-max-no-call-fraction 0.1 \
    }

    output {
        File bam = outputBamBasename + ".bam"
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}