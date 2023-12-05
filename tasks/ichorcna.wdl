version 1.0

import "../structs.wdl"

#ichorCNA related tasks

task SizeselectLessOrEqual150 {
    input {
        File inputBam
        String outputPrefix
        String samtoolsModule = "SAMtools"
        Int? memoryGb = "1"
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
        ml ~{samtoolsModule}

        samtools view -h ~{inputBam} | perl -wlane 'print if(m/^@/ ||(not(m/^@/)&& abs($F[8]) <= 150 ))' | samtools view -S -b -h  -  >  ~{outputPrefix}.bam
        samtools index ~{outputPrefix}.bam
    }
    
    output {
        File bam = outputPrefix + '.bam'
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
    }
}

task HmmcopyReadcounter {
    input {
        File inputBam
        File inputBai
        String outputPrefix
        String hmmcopyutilsModule = "hmmcopy_utils"
        String samtoolsModule = "SAMtools"
        File referencefai

        #String chromosomes = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
        Int windowkilobase = 500
        Int? memoryGb = "1"
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
        Int disk = 1 + ceil(size(inputBam, "G")*2.1)
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
        set -e
        set -o pipefail
        #this part filters out the secondary and supplemental alignments and gives a samtools style index(.bam.bai) cuz the tool wont handle picard style (.bai) ...
        (
            echo "filtering"
            ml ~{samtoolsModule}
            samtools view \
                -h -b \
                --exclude-flags 256 \
                --exclude-flags 2048 \
                ~{inputBam} \
            >  ~{outputPrefix}"_filtered.bam"
            samtools index ~{outputPrefix}"_filtered.bam"
        )
        (
            ml ~{hmmcopyutilsModule}

            #ichorCNA is optimised for human (n chromosomes = 23 plus y chromosome) so this should be fine 
            CHROMOSOMES=$(head -n 24 "~{referencefai}" | \
                cut -f 1 |\
                python3 -c "import sys;chrs = [];[ chrs.append(line.rstrip())  for count,line in enumerate(sys.stdin)]; print(','.join(chrs))")
            readCounter \
            --window ~{windowkilobase}000 \
            --quality 20 \
            --chromosome "$CHROMOSOMES" \
            ~{outputPrefix}"_filtered.bam" > "~{outputPrefix}"".""~{windowkilobase}""kb.wig" 2>/dev/null
        )
    }
    
    output {
        File wiggle = outputPrefix + '.'+ windowkilobase + 'kb.wig'
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
    }
}
#fails integration atm maybe grep reference for main chr matching?
#best to run ichorCNA manually.
#WINDOWKB="500"
#$EBROOTRMINBUNDLEMINICHORCNA/ichorCNA/bin/runIchorCNA.R \
#--WIG  $TUMORWIG \
#--normalPanel sorted_hg38_bams/20230601_normalwigs_${WINDOWKB}kb_pon_median.rds \
#--gcWig /groups/umcg-pmb/tmp01//apps/data/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.${WINDOWKB}000.gc.wig \
#--centromere $EBROOTRMINBUNDLEMINICHORCNA/ichorCNA/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
#--id $(basename $TUMORWIG .wig) \
#--outDir ichorcna_${WINDOWKB}kb/ \
#--mapWig $EBROOTRMINBUNDLEMINICHORCNA/ichorCNA/extdata/map_hg38_${WINDOWKB}kb.wig \
#--genomeBuild "hg38" \
#--genomeStyle "NCBI" \
#--estimateScPrevalence FALSE --scStates "c()" \
#--chrs $(echo "c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX')" | perl -wpe 's/\t/,/g') \
#--chrTrain $(echo "c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr20','chr21','chr22')" | perl -wpe 's/\t/,/g') \
#--chrNormalize $(echo "c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX')" | perl -wpe 's/\t/,/g') \
#--normal "c(0.5,0.75,0.85,0.9,0.95)" &> $PWD/sorted_hg38_bams/20230601_ichorcna_chr19_${WINDOWKB}kb/$(basename $TUMORWIG .wig).ichorcna.run.log