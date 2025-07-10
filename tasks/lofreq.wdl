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
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
    }
    String vcfSuffix =  ".vcf.gz"
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command <<<
        set -o pipefail
        #ordering the files on /tmp
        cat ~{write_lines(select_all(inputBam))} \
            ~{write_lines(inputBamIndex)} | \
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
    }
}
task LoFreqSomatic {
    input {
        File inputNormalBam
        File inputNormalBamIndex
        File inputTumorBam
        File inputTumorBamIndex
        Reference reference
        File targetIntervalList
        String outputVcfBasename
        String lofreqModule = "LoFreq"
        Int? memoryGb = "4"
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
    }
    String vcfSuffix =  ".vcf.gz"
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command <<<
        set -o pipefail
        #ordering the files on /tmp
        cat ~{write_lines(select_all([inputTumorBam,inputNormalBam]))} \
            ~{write_lines([inputNormalBamIndex,inputTumorBamIndex])} | \
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
        
        

        lofreq somatic \
            --ref ~{reference.fasta} \
            --bed ./targets.bed \
            --call-indels \
            --outprefix ~{outputVcfBasename} \
            --normal "$TMPDIR/""$(basename "~{inputBam}")" \
            --tumor "$TMPDIR/""$(basename "~{inputBam}")" 
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
task LoFreqViterbi {
    input {
        File inputBam
        File inputBamIndex
        String outputBasename
        Reference reference
        String lofreqModule = "LoFreq"
        Int? memoryGb = "4"
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command <<<
        set -o pipefail
        #ordering the files on /tmp
        cat ~{write_lines(select_all([inputBam]))} \
            ~{write_lines([inputBamIndex])} | \
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
        
        

        lofreq viterbi \
            --ref ~{reference.fasta} \
            "$TMPDIR/""$(basename "~{inputBam}")"
            --out - | \
        samtools sort - -T $TMPDIR/aln.sorted -o ~{outputBasename}.bam
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
