version 1.0

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
            --output=~{outputBamBasename} \
            --raw-tag=RX \
            --min-map-q=10 \
            --edits=1
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