version 1.0

import "../structs.wdl"

#ichorCNA related tasks

task LoFreqCall {
    input {
        File inputBam
        File inputBamIndex

        Reference reference
        File targetIntervalList
        String outputVcfBasename
        String lofreqModule = "LoFreq"
        Int? memoryGb = "4"
        Int targetScatter = 1
        Int disk = ceil(size([inputBam,inputBamIndex], "M")*1.2)
        Int timeMinutes = 5 + ceil(size(inputBam, "G")) * 120 / targetScatter 
    }
    String vcfSuffix =  ".vcf.gz"
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command <<<
        set -eo pipefail
        #ordering the files on /tmp
        cat ~{write_lines(select_all([inputBam]))} \
            ~{write_lines(select_all([inputBamIndex]))} | \
            (while read FILE; do 
                if [ ! -e "$TMPDIR/""$(basename "$FILE")" ]; then
                    cp "$FILE" "$TMPDIR/"
                    if [[ $FILE =~ \.bam$ ]]; then
                        echo "$TMPDIR/""$(basename "$FILE")">>"./bams_inputs.list" 
                    fi
                else
                    echo "Duplicate file basename spotted $FILE" && exit 1 
                fi
            done )
        ml ~{lofreqModule}
        #on the fly to bed conversion
        grep -v '^@' ~{targetIntervalList} | perl -wlane 'print join("\t",($F[0],$F[1]-1,$F[2],$.));' > targets.bed
        
        

        lofreq call \
            --ref ~{reference.fasta} \
            --bed ./targets.bed \
            --out - \
            "$TMPDIR/""$(basename "~{inputBam}")"| \
            bgzip -c >  \
            ~{outputVcfBasename}~{vcfSuffix} 
    >>>
    
    output {
        File vcf = outputVcfBasename + vcfSuffix
        File vcfIdx = outputVcfBasename + vcfSuffix + ".tbi"
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}
task LoFreqSomatic {
    input {
        File inputNormalBam
        File inputNormalBamIndex
        File inputTumorBam
        File inputTumorBamIndex
        String tumorSampleName
        String normalSampleName
        Reference reference
        File targetIntervalList
        String outputVcfBasename
        String lofreqModule = "LoFreq"
        String picardModule = "picard"
        String pipelineUtilModule= "pipeline-util"
        Int targetScatter = 1
        Int? memoryGb = "4"
        Int timeMinutes = 1 + ceil(size([inputNormalBam,inputTumorBam], "G")) * 120 / targetScatter 
        Int disk = ceil(size([inputNormalBam,inputTumorBam], "M")*1.2)
    }
    String vcfSuffix =  ".vcf.gz"
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command <<<
        set -eo pipefail
        #ordering the files on /tmp
        >&2 echo " ## "$(date)" ## Localising files" 

        cat ~{write_lines(select_all([inputTumorBam,inputNormalBam]))} \
            ~{write_lines(select_all([inputNormalBamIndex,inputTumorBamIndex]))} | \
            (while read FILE; do 
                if [ ! -e "$TMPDIR/""$(basename "$FILE")" ]; then
                    cp "$FILE" "$TMPDIR/"
                    if [[ $FILE =~ \.bam$ ]]; then
                        echo "$TMPDIR/""$(basename "$FILE")">>"./bams_inputs.list" 
                    fi
                else
                    echo "Duplicate file basename spotted $FILE" && exit 1 
                fi
            done )
        
        #on the fly to bed conversion
        grep -v '^@' ~{targetIntervalList} | perl -wlane 'print join("\t",($F[0],$F[1]-1,$F[2],$.));' > targets.bed
        
        (
            ml ~{lofreqModule}
        
            >&2 echo " ## "$(date)" ## Running lofreq somatic" 

            (  
                #if no variants found lofreq will report a non zero exit code.
                lofreq somatic \
                --ref ~{reference.fasta} \
                --bed ./targets.bed \
                --call-indels \
                --outprefix ~{outputVcfBasename} \
                --normal "$TMPDIR/""$(basename "~{inputNormalBam}")" \
                --tumor "$TMPDIR/""$(basename "~{inputTumorBam}")" || \
                (for extension in tumor_stringent.indels.vcf.gz tumor_stringent.snvs.vcf.gz somatic_final.snvs.vcf.gz somatic_final.indels.vcf.gz; do
                    touch  ~{outputVcfBasename}${extension}
                    touch  ~{outputVcfBasename}${extension}.failed
                done )
            )
        )
        #on the fly to bed conversion
        (
            ml ~{pipelineUtilModule}

            >&2 echo " ## "$(date)" ## fixing vcf files" 


            #fixup by moving stuff to the genotype fields and adding general fields (AD/GT)
            for extension in tumor_stringent.indels.vcf.gz tumor_stringent.snvs.vcf.gz somatic_final.snvs.vcf.gz somatic_final.indels.vcf.gz; do
                OUTVCF="~{outputVcfBasename}""$(basename "${extension}" .vcf.gz)"".fixed.vcf.gz"
                if [ $(grep -c '^#' ~{outputVcfBasename}${extension} ) -eq 0 ]; then 
                    cat << EOF | bgzip -c >  ~{outputVcfBasename}$(basename ${extension} .vcf.gz).fixed.vcf.gz
##fileformat=VCFv4.0
##fileDate=20250715
##source=lofreq call -d 101000 -f /groups/umcg-pmb/tmp02/projects/hematopathology/users/umcg-mterpstra/git/Bestie/tests/runs/bamToVariants/cromwell-executions/BamsToVariants/b30f4eea-cde4-444e-93e9-e45dd9ad6f55/call-ScatterAt117_17/shard-0/ScatterAt117_17/034d5aa4-5a46-46ee-81f9-63cb02e75ab8/call-ScatterAt132_21/shard-0/ScatterAt132_21/6b4fb6d9-b76f-4bec-9458-0dfdc3117ac7/call-loFreqSomaticScat/shard-1/inputs/-1706296754/ref.fasta --verbose --no-default-filter -b 1 -l ./targets.bed --call-indels -a 0.010000 -C 7 -s -S lofreq_sample_MySampleS2_norm_MySample_scat0normal_stringent.snvs.vcf.gz,lofreq_sample_MySampleS2_norm_MySample_scat0normal_stringent.indels.vcf.gz -o lofreq_sample_MySampleS2_norm_MySample_scat0tumor_relaxed.vcf.gz /groups/umcg-pmb/tmp02/projects/hematopathology/users/umcg-mterpstra/git/Bestie/tests/runs/bamToVariants/cromwell-executions/BamsToVariants/b30f4eea-cde4-444e-93e9-e45dd9ad6f55/call-ScatterAt117_17/shard-0/ScatterAt117_17/034d5aa4-5a46-46ee-81f9-63cb02e75ab8/call-ScatterAt132_21/shard-0/Sca 
##reference=/groups/umcg-pmb/tmp02/projects/hematopathology/users/umcg-mterpstra/git/Bestie/tests/runs/bamToVariants/cromwell-executions/BamsToVariants/b30f4eea-cde4-444e-93e9-e45dd9ad6f55/call-ScatterAt117_17/shard-0/ScatterAt117_17/034d5aa4-5a46-46ee-81f9-63cb02e75ab8/call-ScatterAt132_21/shard-0/ScatterAt132_21/6b4fb6d9-b76f-4bec-9458-0dfdc3117ac7/call-loFreqSomaticScat/shard-1/inputs/-1706296754/ref.fasta
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
##INFO=<ID=SB,Number=1,Type=Integer,Description="Phred-scaled strand bias at this position">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=CONSVAR,Number=0,Type=Flag,Description="Indicates that the variant is a consensus variant (as opposed to a low frequency variant).">
##INFO=<ID=HRUN,Number=1,Type=Integer,Description="Homopolymer length to the right of report indel position">
##FILTER=<ID=min_dp_7,Description="Minimum Coverage 7">
##FILTER=<ID=max_dp_100000,Description="Maximum Coverage 100000">
##FILTER=<ID=sb_fdr,Description="Strand-Bias Multiple Testing Correction: fdr corr. pvalue > 0.001000">
##FILTER=<ID=indelqual_bonf,Description="Indel Quality Multiple Testing Correction: bonf corr. pvalue < 0.010000">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	~{tumorSampleName}	~{normalSampleName}
EOF
                else
                    #This converts the lofreq output to a cleaner formatted vcf file with the tumor sample speficic data in the info fields
                    #This also adds in the normalsample as a filler for completeness.  
                    perl -wpe 's/^##fileformat=VCFv4\.0$/##fileformat=VCFv4.2/;
                        if($_ =~ /^#CHROM/){
                            print "##INFO=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
                            print "##INFO=<ID=AD,Number=2,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n";
                            print "##INFO=<ID=UNIQ,Number=0,Type=Flag,Description=\"Unique, i.e. not detectable in paired sample\">\n##INFO=<ID=UQ,Number=1,Type=Integer,Description=\"Phred-scaled uniq score at this position\">\n";
                            print "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic event\">\n";
                        }
                        if(m/DP4=(\d+(,\d+){3})(;|\n)/){
                            my @dp4=split(",",$1);
                            $_=substr($_,0,-1).";AD=".($dp4[0]+$dp4[1]).",".($dp4[2]+$dp4[3])."\n";
                            if($dp4[0]+$dp4[1]> 0){
                                $_=substr($_,0,-1).";GT=0/1\n";
                            }else{
                                $_=substr($_,0,-1).";GT=1/1\n";
                            }
                        };' <(bgzip -dc ~{outputVcfBasename}${extension} ) > ~{outputVcfBasename}$(basename ${extension} .vcf.gz).tmp.fixed.vcf
                    InfoFieldsToGenotypeFields.pl  \
                    -f 'GT,AD,DP,DP4,AF,SB' \
                    -g "~{tumorSampleName}" \
                    -i ~{outputVcfBasename}$(basename ${extension} .vcf.gz).tmp.fixed.vcf | \
                    perl -wpe 's/~{tumorSampleName}/~{tumorSampleName}\t~{normalSampleName}/ if(m/^#CHROM/); s/\n/\t.\n/ if(m/^[^#]/)' | \
                    bgzip -c > ~{outputVcfBasename}$(basename ${extension} .vcf.gz).fixed.vcf.gz;
                fi
            done
        )
        >&2 echo " ## "$(date)" ## picard sorting vcf files" 

        #Merge tumor calls to one VCF 
        (
            ml ~{picardModule}
            java -jar $EBROOTPICARD/picard.jar SortVcf \
                SD=~{reference.dict} \
                I=~{outputVcfBasename}tumor_stringent.indels.fixed.vcf.gz I=~{outputVcfBasename}tumor_stringent.snvs.fixed.vcf.gz \
                O=~{outputVcfBasename}tumor_stringent.fixed.vcf.gz
        )
        >&2 echo " ## "$(date)" ## done" 

    >>>
    
    output {
        File vcf = outputVcfBasename + 'tumor_stringent.fixed' + vcfSuffix
        File vcfIdx = vcf + ".tbi"
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}
task LoFreqViterbi {
    input {
        IndexedFile inputBamIndexed
        String outputBasename
        Reference reference
        String lofreqModule = "LoFreq"
        String picardModule = "picard"
        Int? memoryGb = "8"
        Int timeMinutes = 1 + ceil(size(inputBamIndexed.file, "G")) * 120
        Int disk = ceil(size(inputBamIndexed.file, "M")*1.2)  
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command <<<
        set -eo pipefail
        #ordering the files on /tmp
        cat ~{write_lines(select_all([inputBamIndexed.file]))} \
            ~{write_lines([inputBamIndexed.index])} | \
            (while read FILE; do 
                if [ ! -e "$TMPDIR/""$(basename "$FILE")" ]; then
                    cp "$(realpath "$FILE")" "$TMPDIR/"
                    if [[ $FILE =~ \.bam$ ]]; then
                        echo "$TMPDIR/""$(basename "$FILE")">>"./bams_inputs.list" 
                    fi
                else
                    echo "Duplicate file basename spotted $FILE" && exit 1 
                fi
            done )
        ml ~{lofreqModule}
        if [ ! -e ~{inputBamIndexed.file}.bai ]; then 
            cp -v "$TMPDIR/""$(basename "~{inputBamIndexed.index}")" "$TMPDIR/""$(basename "~{inputBamIndexed.file}.bai")"
        fi
        lofreq viterbi \
            --ref "~{reference.fasta}" \
            "$TMPDIR/""$(basename "~{inputBamIndexed.file}")" \
            --out - | \
        
        (
            ml ~{picardModule} && java -jar $EBROOTPICARD/picard.jar \
            FixMateInformation \
            SORT_ORDER=coordinate \
            I=/dev/stdin O=/dev/stdout \
        ) | (
            ml ~{picardModule} && java -jar $EBROOTPICARD/picard.jar \
            SetNmMdAndUqTags \
            REFERENCE_SEQUENCE="~{reference.fasta}" \
            I=/dev/stdin O=~{outputBasename}.bam CREATE_INDEX=true
        )
        cp ~{outputBasename}.bai ~{outputBasename}.bam.bai  

    >>>
    
    output {
        File bam = outputBasename + '.bam'
        File bai = outputBasename + '.bam.bai'
        IndexedFile bamOut = { 
          "file" : outputBasename + '.bam',
          "index" : outputBasename + ".bam.bai"
        }
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}
task LoFreqViterbiCached {
    input {
        IndexedFile inputBamIndexed
        String outputBasename
        Reference reference
        String lofreqModule = "LoFreq"
        String picardModule = "picard"
        Int? memoryGb = "8"
        Int timeMinutes = 1 + ceil(size(inputBamIndexed.file, "G")) * 120
        Int disk = ceil(size(inputBamIndexed.file, "M")*1.2)  
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command <<<
        set -xeo pipefail
        if [ -e "/groups/umcg-pmb/tmp02/projects/hematopathology/users/umcg-mterpstra/projects/20250724_braf_S6_variant_analysis/lofreq-cache/""$(basename "~{outputBasename}.bam")" ]; then 
            
            >&2 echo " ## "$(date)" ## Found cached result copying and skipping" 


            cp "/groups/umcg-pmb/tmp02/projects/hematopathology/users/umcg-mterpstra/projects/20250724_braf_S6_variant_analysis/lofreq-cache/""$(basename "~{outputBasename}.bam")" ~{outputBasename}.bam
            cp "/groups/umcg-pmb/tmp02/projects/hematopathology/users/umcg-mterpstra/projects/20250724_braf_S6_variant_analysis/lofreq-cache/""$(basename "~{outputBasename}.bam")".bai ~{outputBasename}.bam.bai
            exit 0
        fi
        #ordering the files on /tmp

        >&2 echo " ## "$(date)" ## Localising files" 

        cat ~{write_lines(select_all([inputBamIndexed.file]))} \
            ~{write_lines([inputBamIndexed.index])} | \
            (while read FILE; do 
                if [ ! -e "$TMPDIR/""$(basename "$FILE")" ]; then
                    cp "$FILE" "$TMPDIR/"
                    if [[ $FILE =~ \.bam$ ]]; then
                        echo "processing bam file $TMPDIR/""$(basename "$FILE")"
                        (ml SAMtools && 
                            samtools view -h $FILE | head -n 200000 > test.sam
                            samtools view -Shb test.sam >  "$TMPDIR/""$(basename "$FILE")"
                            samtools index "$TMPDIR/""$(basename "$FILE")"
                        )
                        echo "$TMPDIR/""$(basename "$FILE")">>"./bams_inputs.list" 
                    elif [[ $FILE =~ \.bai$ ]]; then
                        echo "Skipped bai file $TMPDIR/""$(basename "$FILE")"
                    else
                        cp "$FILE" "$TMPDIR/"
                    fi
                else
                    echo "Duplicate file basename spotted $FILE" && exit 1 
                fi
            done )
        >&2 echo " ## "$(date)" ## Loading module" 

        ml ~{lofreqModule}
        if [ ! -e ~{inputBamIndexed.file}.bai ]; then 
            cp -v "$TMPDIR/""$(basename "~{inputBamIndexed.index}")" "$TMPDIR/""$(basename "~{inputBamIndexed.file}.bai")"
        fi
        >&2 echo " ## "$(date)" ## Running lofreq viterbi and fixing up output" 

        lofreq viterbi \
            --ref "~{reference.fasta}" \
            "$TMPDIR/""$(basename "~{inputBamIndexed.file}")" \
            --out - | \
        
        (
            ml ~{picardModule} && java -jar $EBROOTPICARD/picard.jar \
            FixMateInformation \
            SORT_ORDER=coordinate \
            I=/dev/stdin O=/dev/stdout \
        ) | (
            ml ~{picardModule} && java -jar $EBROOTPICARD/picard.jar \
            SetNmMdAndUqTags \
            REFERENCE_SEQUENCE="~{reference.fasta}" \
            I=/dev/stdin O=~{outputBasename}.bam CREATE_INDEX=true
        )
        cp ~{outputBasename}.bai ~{outputBasename}.bam.bai

        >&2 echo " ## "$(date)" ## Lofreq viterbi done " 


    >>>
    
    output {
        File bam = outputBasename + '.bam'
        File bai = outputBasename + '.bam.bai'
        IndexedFile bamOut = { 
          "file" : outputBasename + '.bam',
          "index" : outputBasename + ".bam.bai"
        }
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}
