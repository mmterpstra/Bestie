version 1.0

task fqstats {
    input {
        File fastq_input
        Int mem_gb
	String fastqc_mod
    }

    command {
        module load ${fastqc_mod} && fastqc -o ./ ${fastq_input} 
    }

    output {
        File outHtml = basename(basename(fastq_input, ".fastq.gz"), ".fq.gz") + "_fastqc.html"
        File outZip = basename(basename(fastq_input, ".fastq.gz"), ".fq.gz") + "_fastqc.zip"
    }

    runtime {
        memory: mem_gb + "GB"
    }
}

workflow fastqstatsWorkflow {
    input {
        File fastq_input
        Int mem_gb
	String fastqc_mod
    }
    
    call fqstats { input: fastq_input=fastq_input, mem_gb=mem_gb, fastqc_mod=fastqc_mod }

    meta {
        author: "MMTerpstra"
        description: "## Fastqc stats \n This is the fastqc stats workflow."
    }
}

