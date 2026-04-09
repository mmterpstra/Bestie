version 1.0

import "../structs.wdl" as structs

task RemoveNumis {
    input {
        File inputBam
        String outputBasename
        Int memoryGb = "1"
        String samtoolsModule = "SAMtools"
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
        #Int disk = 1 + ceil(size(inputBam, "G")) * 1024
    }
    command <<<
        set -e -o pipefail
        module load ~{samtoolsModule} && \
        (
            samtools view -h ~{inputBam} | \
            perl -wne 'BEGIN{$|=1};print if(m/^\@/ or not(m/N[ATCG]*\-/ or m/\-[ATCG]*N/))' | \
            samtools view -h -b -1 
        )> ~{outputBasename}.bam
 

    >>>

    output {
        File bam = outputBasename + ".bam"

    }
}
#wip
task Stats {
    input {
        File inputBam
        File inputBai
        Reference reference
        String outputBasename
        Int memoryGb = "1"
        String samtoolsModule = "SAMtools"
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
        #Int disk = 1 + ceil(size(inputBam, "G")) * 1024
    }
    command <<<
        set -e -o pipefail
        module load ~{samtoolsModule} && \
        (
            samtools stats --ref-seq "~{reference.fasta}" "~{inputBam}" > "~{outputBasename}.samstats"
        )
 

    >>>

    output {
        File stats = outputBasename + ".samstats"

    }
}
task XYIdxStats {
    input {
        File inputBam
        File inputBai
        String outputBasename
        Int memoryGb = "1"
        String samtoolsModule = "SAMtools"
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
        #Int disk = 1 + ceil(size(inputBam, "G")) * 1024
    }
    command <<<
        set -e -o pipefail
        module load ~{samtoolsModule} && \
        (
            samtools idxstats  "~{inputBam}" > "~{outputBasename}.tmp"
            head -n 24 "~{outputBasename}.tmp" | tail -n 2 > "~{outputBasename}.samidxstats"
            rm -v "~{outputBasename}.tmp"
        )
    >>>

    output {
        File idxstats = outputBasename + ".samidxstats"

    }
}


