#test/data/
This dir contains test data either validation or for input needed to run the pipeline.

#fastq.gz inputs

Creating fastq inputs and fasta reference for analysis

```
(export REFFA=/groups/umcg-pmb/tmp01/apps/data/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
    ml BEDTools/2.27.1-foss-2018b wgsim/a12da33-foss-2018b
    bedtools getfasta -fi $REFFA \
        -bed <( echo -e "chr8\t127736594\t127740958\tNMYC") \
        -fo /dev/stdout > tests/data/ref/ref.fasta
    wgsim tests/data/ref/ref.fasta tests/data/raw/fastq/reads_R1.fq tests/data/raw/fastq/reads_R2.fq \
        &>tests/data/raw/fastq/reads.wgsim.log
    gzip  tests/data/raw/fastq/reads_R1.fq tests/data/raw/fastq/reads_R2.fq)
(mv tests/data/raw/fastq/reads_R1.fq.gz{,.bak}
    mv tests/data/raw/fastq/reads_R2.fq.gz{,.bak};
    zcat tests/data/raw/fastq/reads_R1.fq.gz.bak  | perl -wpe 'if($.%4==0){chomp;print "=" x length($_)."\n"; $_=""}' | gzip -c > tests/data/raw/fastq/reads_R1.fq.gz;
    zcat tests/data/raw/fastq/reads_R2.fq.gz.bak  | perl -wpe 'if($.%4==0){chomp;print "=" x length($_)."\n"; $_=""}' | gzip -c > tests/data/raw/fastq/reads_R2.fq.gz)
(ml SAMtools/1.16.1-GCCcore-11.3.0; samtools faidx tests/data/ref/ref.fasta)
(ml GATK/4.2.4.1-Java-8-LTS; gatk CreateSequenceDictionary --REFERENCE tests/data/ref/ref.fasta )
```
```
#bwa index creation
(ml BWA/0.7.17-GCCcore-11.3.0 &&  bwa index tests/data/ref/ref.fasta)
#unmapped sam/bam
(ml picard/2.26.10-Java-8-LTS && java -Xmx8g -jar $EBROOTPICARD/picard.jar FastqToSam FASTQ=tests/data/raw/fastq/reads_R1.fq.gz FASTQ2=tests/data/raw/fastq/reads_R2.fq.gz SAMPLE_NAME="SAMPLE" PLATFORM="ILLUMINA" OUTPUT=tests/data/bam/unaligned_read_pairs.bam)
#aligned bam unbam>interleaved.fq>bwa.sam>merged.bam
(ml SAMtools && samtools faidx tests/data/ref/ref.fasta)
(ml picard/2.26.10-Java-8-LTS &&  java -Xms1000m -Xmx1000m  -jar $EBROOTPICARD/picard.jar \
	CreateSequenceDictionary REFERENCE=tests/data/ref/ref.fasta
    
    java -Xms1000m -Xmx1000m  -jar $EBROOTPICARD/picard.jar \
	SamToFastq \
	INPUT=tests/data/bam/unaligned_read_pairs.bam \
	FASTQ=/dev/stdout \
	INTERLEAVE=true \
	NON_PF=true | \
	(ml purge && ml BWA/0.7.17-GCCcore-11.3.0 &&  bwa mem -K 100000000 -p -v 3 -t 16 -Y tests/data/ref/ref.fasta /dev/stdin - 2> >(tee tests/data/ref/ref.fasta.bwa.stderr.log >&2)) | \
	java -Dsamjdk.compression_level=1 -Xms1000m -Xmx1000m -jar $EBROOTPICARD/picard.jar \
    MergeBamAlignment \
        VALIDATION_STRINGENCY=SILENT \
        EXPECTED_ORIENTATIONS=FR \
        ATTRIBUTES_TO_RETAIN=X0 \
        ATTRIBUTES_TO_REMOVE=NM \
        ATTRIBUTES_TO_REMOVE=MD \
        ALIGNED_BAM=/dev/stdin \
        UNMAPPED_BAM=tests/data/bam/unaligned_read_pairs.bam \
        OUTPUT=tests/data/bam/aligned_read_pairs.bam \
        REFERENCE_SEQUENCE=tests/data/ref/ref.fasta \
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
        ADD_PG_TAG_TO_READS=false
)

```