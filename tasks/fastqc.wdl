version 1.0

task FastQC {
    input {
        File inputFastq
        Int? memoryGb = "1"
	    String fastqcModule
        Int timeMinutes = 1 + ceil(size(inputFastq, "G")) * 20
    }

    command {
        set -e
        module load ~{fastqcModule} && fastqc -o ./ ~{inputFastq}
    }

    output {
        String fastqcDir =  basename(inputFastq, ".fastq.gz")
        #ile out_html1 = basename(inputFastq1, ".fastq.gz") + "_fastqc.html"
        File outZip = basename(inputFastq, ".fastq.gz") + "_fastqc.zip"
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
    }
}
task FastQCSample {
    input {
        Array[File] inputFastqGzs
        String outputPrefix
        Int? memoryGb = "1"
	    String fastqcModule
        Int timeMinutes = 1 + ceil(size(inputFastqGzs, "G")) * 20
    }

    command <<<
        set -e -o pipefail
        module load ~{fastqcModule}
        mkfifo ~{outputPrefix}.fq
        zcat ~{sep=" " inputFastqGzs} > ~{outputPrefix}.fq &
        fastqc -o ./ ~{outputPrefix}.fq
    >>>

    output {
        #File out_html1 = basename(inputFastq1, ".fastq.gz") + "_fastqc.html"
        File outZip = outputPrefix + "_fastqc.zip"
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
    }
}
task FastQCSamplePaired {
    input {
        Array[File] inputFastqGzs
        Array[File?] inputFastq2Gzs
        String outputPrefix
        Int? memoryGb = "1"
	    String fastqcModule
        Int timeMinutes = 1 + ceil(size(inputFastqGzs, "G")) * 20
    }

    command <<<
        set -e -o pipefail
        module load ~{fastqcModule}
        mkfifo ~{outputPrefix}_R1.fq
        zcat ~{sep=" " inputFastqGzs} > ~{outputPrefix}_R1.fq &
        if [[ -f ~{sep=" " inputFastq2Gzs[0]} ]] ; then
            mkfifo ~{outputPrefix}_R2.fq
            zcat ~{sep=" " inputFastq2Gzs} > ~{outputPrefix}_R2.fq &
            fastqc -o ./ ~{outputPrefix}.fq ~{outputPrefix}_R2.fq
        else
            fastqc -o ./ ~{outputPrefix}.fq
        fi
    >>>

    output {
        #File out_html1 = basename(inputFastq1, ".fastq.gz") + "_fastqc.html"
        File outZip = outputPrefix + "_fastqc.zip"
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
    }
}
task FastQCPaired {
    input {
        File inputFastq
        File? inputFastq2
        Int? memoryGb = "1"
	    String fastqcModule
        Int timeMinutes = 1 + ceil(size(inputFastq, "G")) * 20
        String outputFastqcBasename = "sample_fastqc"
    }

    command {
        set -e
        module load ~{fastqcModule}
        mkdir -p ~{outputFastqcBasename}/
        fastqc -o ~{outputFastqcBasename} ~{inputFastq} ~{inputFastq2}
        zip -ru ~{outputFastqcBasename}.zip ~{outputFastqcBasename}
    }

    output {
        #String fastqcDir =  basename(inputFastq, ".fastq.gz")
        #ile out_html1 = basename(inputFastq1, ".fastq.gz") + "_fastqc.html"
        File zip = outputFastqcBasename + ".zip"
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
    }
}