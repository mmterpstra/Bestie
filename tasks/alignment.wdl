version 1.0

#ill put the minimalist unaligned.bam to algned bam here for optimalisation
import "../structs.wdl" as structs

task bwaAlignBam {
    input {
        File inputUnalignedBam
        #note on memory usage
        #mainly ~8g for mergeBamAlignment memory, ~8g for bwa/bamtofastq and optionally ~8 for sorting/extra umiconsensus read data overhead.
        #This is for 150 bp reads for longer reads maybe more?
        Int memoryGb = if coordinateSort then 16+8 else 16
        #optional filling in only works partailly idk why
        Int javaMemoryMb = if coordinateSort then ceil((16+8-8)*1024*0.95) else ceil((16-8)*1024*0.95)
	    String bwaModule = "BWA/0.7.17-GCCcore-11.3.0"
        String picardModule = "picard/2.26.10-Java-8-LTS"
        BwaIndex referenceBwaIndex
        Reference reference
        Boolean umiTags = false
        String outputBamBasename
        Int threads = 8
        Int disk = ceil(size(inputUnalignedBam, "M")*2.1)
        Int timeMinutes = 1 + ceil(size(inputUnalignedBam, "G")) * 120 + 20
        Boolean coordinateSort = false
    }
    #Int javaMemoryMb =  if coordinateSort then ceil((memoryGb+8-8) * 1024) else ceil(memoryGb-8 * 1024)
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
            -p -v 3 -t ~{threads - 1} -Y \
            ~{referenceBwaIndex.fastaFile} \
            /dev/stdin - 2> >(tee ./bwa.stderr.log >&2)) | \
        (ml load ~{picardModule} && \
            java -Dsamjdk.compression_level=1 -Xms1000m -Xmx~{javaMemoryMb}m -jar $EBROOTPICARD/picard.jar \
            MergeBamAlignment \
            VALIDATION_STRINGENCY=SILENT \
            EXPECTED_ORIENTATIONS=FR \
            ATTRIBUTES_TO_RETAIN=X0 \
            ATTRIBUTES_TO_REMOVE=NM \
            ATTRIBUTES_TO_REMOVE=MD \
            ALIGNED_BAM=/dev/stdin \
            UNMAPPED_BAM=~{inputUnalignedBam} \
            OUTPUT=~{outputBamBasename}.bam \
            REFERENCE_SEQUENCE=~{reference.fasta} \
            SORT_ORDER=~{if coordinateSort then "\"coordinate\"" else "\"unsorted\""} \
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
            PROGRAM_GROUP_COMMAND_LINE="bwa mem -K 100000000 -p -v 3 -t ~{threads - 1} -Y" \
            PROGRAM_GROUP_NAME="bwamem" \
            UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
            ALIGNER_PROPER_PAIR_FLAGS=true \
            ADD_PG_TAG_TO_READS=false \
            ~{if umiTags then "ATTRIBUTES_TO_RETAIN=ZS" else ""} \
            ~{if umiTags then "ATTRIBUTES_TO_RETAIN=ZI" else ""} \
            ~{if umiTags then "ATTRIBUTES_TO_RETAIN=ZM" else ""} \
            ~{if umiTags then "ATTRIBUTES_TO_RETAIN=ZC" else ""} \
            ~{if umiTags then "ATTRIBUTES_TO_RETAIN=ZN" else ""} \
            ~{if umiTags then "ATTRIBUTES_TO_RETAIN=ad" else ""} \
            ~{if umiTags then "ATTRIBUTES_TO_RETAIN=bd" else ""} \
            ~{if umiTags then "ATTRIBUTES_TO_RETAIN=cd" else ""} \
            ~{if umiTags then "ATTRIBUTES_TO_RETAIN=ae" else ""} \
            ~{if umiTags then "ATTRIBUTES_TO_RETAIN=be" else ""} \
            ~{if umiTags then "ATTRIBUTES_TO_RETAIN=ce" else ""} \
            ~{if coordinateSort then "CREATE_INDEX=true" else ""} 
            )
    >>>

    output {
        File bam = outputBamBasename + ".bam"
        File? bai = outputBamBasename + ".bai"
    }

    runtime {
        memory: ceil(memoryGb * 1024)
        timeMinutes: timeMinutes
        cpus: threads
        disk: disk 
    }
}