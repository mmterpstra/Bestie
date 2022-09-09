version 1.0

task bamstats {
    input {
        File bam_input
        Int mem_gb
    }

    command {
        echo ${mem_gb} ${bam_input} > bamstats_report.txt
    }

    output {
        File bamstats_report = "bamstats_report.txt"
    }

    runtime {
        memory: mem_gb + "GB"
    }
}

workflow bamstatsWorkflow {
    input {
        File bam_input
        Int mem_gb

    }
    
    call bamstats { input: bam_input=bam_input, mem_gb=mem_gb }

    meta {
        author: "Andrew Duncan"
        email: "andrew@foobar.com"
        description: "## Bamstats \n This is the Bamstats workflow.\n\n Adding documentation improves clarity."
    }
}

