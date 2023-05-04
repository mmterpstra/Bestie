version 1.0

#ill put the minimalist unaligned.bam to algned bam here for optimalisation
import "../structs.wdl" as structs

task bwaAlignBam {
    input {
        File inputUnalignedBam
        Int? memoryGb = "12"
	    String bwaModule = "BWA/0.7.17-GCCcore-11.3.0"
        String picardModule = "picard/2.26.10-Java-8-LTS"
        BwaIndex referenceBwaIndex
        Reference reference
        String outputBam
        Int timeMinutes = 1 + ceil(size(inputUnalignedBam, "G")) * 40
    }
    command <<<
        set -o pipefail

        (ml load ~{picardModule} && \
            java -Xms1000m -Xmx1000m  -jar $EBROOTPICARD/picard.jar \
            SamToFastq \
            INPUT=~{inputUnalignedBam} \
            FASTQ=/dev/stdout \
            INTERLEAVE=true \
            NON_PF=true) | \
        (ml load ~{bwaModule} && \
            bwa mem -K 100000000 \
            -p -v 3 -t 16 -Y \
            ~{referenceBwaIndex.fastaFile} \
            /dev/stdin - 2> >(tee ./bwa.stderr.log >&2)) | \
        (ml load ~{picardModule} && \
            java -Dsamjdk.compression_level=1 -Xms1000m -Xmx1000m -jar $EBROOTPICARD/picard.jar \
            MergeBamAlignment \
            VALIDATION_STRINGENCY=SILENT \
            EXPECTED_ORIENTATIONS=FR \
            ATTRIBUTES_TO_RETAIN=X0 \
            ATTRIBUTES_TO_REMOVE=NM \
            ATTRIBUTES_TO_REMOVE=MD \
            ALIGNED_BAM=/dev/stdin \
            UNMAPPED_BAM=~{inputUnalignedBam} \
            OUTPUT=~{outputBam} \
            REFERENCE_SEQUENCE=~{reference.fasta} \
            SORT_ORDER="unsorted" \
            IS_BISULFITE_SEQUENCE=false \
            ALIGNED_READS_ONLY=false \
            CLIP_ADAPTERS=false \
            CLIP_OVERLAPPING_READS=true \
            MAX_RECORDS_IN_RAM=2000000 \
            ADD_MATE_CIGAR=true \
            MAX_INSERTIONS_OR_DELETIONS=-1 \
            PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
            PROGRAM_RECORD_ID="bwamem" \
            PROGRAM_GROUP_VERSION="0.7.17-r1188" \
            PROGRAM_GROUP_COMMAND_LINE="bwa mem -K 100000000 -p -v 3 -t 16 -Y" \
            PROGRAM_GROUP_NAME="bwamem" \
            UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
            ALIGNER_PROPER_PAIR_FLAGS=true \
            ADD_PG_TAG_TO_READS=false)
    >>>

    output {
        File alignedBam = outputBam
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
    }
}