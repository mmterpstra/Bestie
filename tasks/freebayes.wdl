version 1.0

import "../structs.wdl"

#ichorCNA related tasks

task Freebayes {
    input {
        Array [File] inputBams
        Array [File] inputBamIndexes

        Reference reference
        File targetIntervalList
        String outputVcfBasename
        String freebayesModule = "freebayes"
        Int targetScatter = 1
        Int? memoryGb = "18"
        Int timeMinutes = 1 + ceil(size(inputBams, "G")) * 120 / targetScatter 
        
    }
    String vcfSuffix =  ".vcf.gz"
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command <<<
        set -o pipefail
        #ordering the files on /tmp
        cat ~{write_lines(select_all(inputBams))} \
            ~{write_lines(inputBamIndexes)} | \
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
        ml ~{freebayesModule}
        #on the fly to bed conversion
        grep -v '^@' ~{targetIntervalList} | perl -wlane 'print join("\t",($F[0],$F[1]-1,$F[2],$.));' > targets.bed
        
        
        freebayes \
            --fasta-reference ~{reference.fasta} \
            --targets ./targets.bed \
            --bam-list ./bams_inputs.list| \
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
    }
}

task FreebayesSomatic {
    input {
        Array [File] inputBams
        Array [File] inputBamIndexes
        Reference reference
        File targetIntervalList
        String outputVcfBasename
        String freebayesModule = "freebayes"
        Int? memoryGb = "18"
        Int targetScatter = 1
        Int timeMinutes = 1 + ceil(size(inputBams, "G")) * 400 /targetScatter
        Int disk = ceil(size(inputBams, "M")*1.2)
    }
    String vcfSuffix =  ".vcf.gz"

    command <<<
        set -e 
        set -o pipefail
        #idk why but freebayes randomises the vcf output file sample columns.  
        >&2 echo " ## "$(date)" ## Localising files" 
        cat ~{write_lines(inputBams)} \
            ~{write_lines(inputBamIndexes)} | sort | \
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
        >&2 echo " ## "$(date)" ## Loading modules" 
        ml ~{freebayesModule}

        >&2 echo " ## "$(date)" ## Converting interval_list to bed" 
        #on the fly to bed conversion
        grep -v '^@' ~{targetIntervalList} | perl -wlane 'print join("\t",($F[0],$F[1]-1,$F[2],$.));' > targets.bed
        
        #this limits the output on wes/wgs
        minAlternateFraction="0.03"
        if [ $(perl -wane 'use POSIX;BEGIN{our $sum;} next if(m/^#|^track/);$sum += $F[2] - $F[1];END{print (floor($sum / 10000000 *20)."\n");}' ./targets.bed ) -ge 1 ];then 
            minAlternateFraction="0.001"
        else 
            minAlternateFraction="0.03"
        fi

        #idk why but depending on amount of intervals/reference retrieved 
        #from remote disk performance may suffer > 10x slowdowns
        
        >&2 echo " ## "$(date)" ## Running freebayes" 
        freebayes \
            --fasta-reference ~{reference.fasta} \
            --allele-balance-priors-off \
            --binomial-obs-priors-off \
            --hwe-priors-off \
            --min-alternate-fraction $minAlternateFraction \
            --min-mapping-quality 20 \
            --max-complex-gap 20 \
            --genotype-qualities \
            --report-genotype-likelihood-max \
            --min-alternate-count 3 \
            --pooled-continuous \
            --strict-vcf \
            --vcf /dev/stdout \
            --targets ./targets.bed \
            --bam-list ./bams_inputs.list \
            --use-best-n-alleles 5 \
            --limit-coverage 10000 \
            --haplotype-length -1 | \
        perl -wpe 's/\t\.:\.(:\.)+/\t./ if m/\t\.:\.(:\.)+[\t\n]/' | \
        tee >( perl -wne 'if($.%20 == 0){
            print "## info ## ".scalar(localtime)." ## ".$_ ;}' \
                >> /dev/stderr ) | \
        bgzip -c >  \
            ~{outputVcfBasename}~{vcfSuffix}
        
        #this should store the index appending .tbi to the filename 

        >&2 echo " ## "$(date)" ## Indexing output"

        tabix -p vcf ~{outputVcfBasename}~{vcfSuffix}

        >&2 echo " ## "$(date)" ## Done"
    >>>
    
    output {
        File vcf = outputVcfBasename + vcfSuffix
        File vcfIdx = outputVcfBasename + vcfSuffix + ".tbi"
        IndexedFile vcfOut = { 
          "file" : vcf,
          "index" : vcfIdx
        }
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}


task FreebayesRecall {
    #freebayes -f ref.fa -@ in.vcf.gz aln.bam >var.vcf
    input {
        Array [File] inputBams
        Array [File] inputBamIndexes
        Reference reference
        IndexedFile inputVariants
        String outputVcfBasename
        String freebayesModule = "freebayes"
        Int? memoryGb = "18"
        Int targetScatter = 1
        Int timeMinutes = 5 + ceil(size(inputBams, "G")) * 400 / targetScatter
        Int disk = ceil(size(inputBams, "M")*1.2)
        
    }
    String vcfSuffix =  ".vcf.gz"

    command <<<
        set -e 
        set -o pipefail

        >&2 echo " ## "$(date)" ## Localising files" 
        cat ~{write_lines(inputBams)} \
            ~{write_lines(inputBamIndexes)} | \
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
        >&2 echo " ## "$(date)" ## Loading modules" 
        ml ~{freebayesModule}
        
        >&2 echo " ## "$(date)" ## Running freebayes" 
        freebayes \
            --fasta-reference ~{reference.fasta} \
            --variant-input ~{inputVariants.file} \
            --only-use-input-alleles \
            --haplotype-length 0 \
            --min-alternate-count 1 \
            --min-alternate-fraction 0 \
            --no-population-priors \
            --allele-balance-priors-off \
            --report-monomorphic \
            --vcf /dev/stdout \
            --bam-list ./bams_inputs.list \
            --limit-coverage 20000 | \
        perl -wpe 's/\t\.:\.(:\.)+/\t./ if m/\t\.:\.(:\.)+[\t\n]/' | \
        tee >( perl -wne 'if($.%20 == 0){
            print "## info ## ".scalar(localtime)." ##".$_ ;}' \
                >> /dev/stderr ) | \
        bgzip -c >  \
            ~{outputVcfBasename}~{vcfSuffix}
        
        #this should store the index appending .tbi to the filename 

        >&2 echo " ## "$(date)" ## Indexing output"

        tabix -p vcf ~{outputVcfBasename}~{vcfSuffix}
        >&2 echo " ## "$(date)" ## Validating output call lengths"
        if [ "$(grep -c -v '^#' ~{inputVariants.file})" -ne "$(grep -c -v '^#' ~{outputVcfBasename}~{vcfSuffix})" ] ; then
            echo "Input and output variant call counts aren't equal. Check '~{inputVariants.file}' and '~{outputVcfBasename}~{vcfSuffix}'." &&  exit 1 
        fi
        >&2 echo " ## "$(date)" ## Done"
    >>>
    
    output {
        File vcf = outputVcfBasename + vcfSuffix
        File vcfIdx = outputVcfBasename + vcfSuffix + ".tbi"
        IndexedFile vcfOut = { 
          "file" : outputVcfBasename + vcfSuffix,
          "index" : outputVcfBasename + vcfSuffix + ".tbi"
        }
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}