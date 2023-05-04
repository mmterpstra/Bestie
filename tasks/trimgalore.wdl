version 1.0

task TrimGalore {
    input {
        File inputFastq1
        File? inputFastq2
        String outputFastq1
        String? outputFastq2
        Int? memoryGb = 1
        Int timeMinutes = 1 + ceil(size(inputFastq1, "G")) * 50
        #File? fastq_input_umi
        #String samplename
        #String identifier
	    String trimgaloreModule
    }
    command {
        set -e 
        module load ${trimgaloreModule} && \
        if [ "${inputFastq2}x" == "x" ];then \
            trim_galore "${inputFastq1}" \
                --output_dir "$(basename $(basename ${outputFastq1} .fastq.gz) .fq.gz)"
                ln -sf $PWD/$(basename $(basename ${outputFastq1} .fastq.gz) .fq.gz)/$(basename $(basename ${inputFastq1} .fastq.gz) .fq.gz)_trimmed.fq.gz ${outputFastq1}
                ln -sf $PWD/$(basename $(basename ${outputFastq1} .fastq.gz) .fq.gz)/$(basename ${inputFastq1} )"_trimming_report.txt" ${outputFastq1}"_trimming_report.txt"
        else \
            trim_galore --paired "${inputFastq1}" "${inputFastq2}" \
              --output_dir "$(basename $(basename ${outputFastq1} .fastq.gz) .fq.gz)" 
            ln -sf $PWD/$(basename $(basename ${outputFastq1} .fastq.gz) .fq.gz)/$(basename $(basename ${inputFastq1} .fastq.gz) .fq.gz)_val_1.fq.gz ${outputFastq1}
            ln -sf $PWD/$(basename $(basename ${outputFastq1} .fastq.gz) .fq.gz)/$(basename ${inputFastq1} )"_trimming_report.txt" ${outputFastq1}"_trimming_report.txt"
            ln -sf $PWD/$(basename $(basename ${outputFastq1} .fastq.gz) .fq.gz)/$(basename $(basename ${inputFastq2} .fastq.gz) .fq.gz)_val_2.fq.gz ${outputFastq2}
            ln -sf $PWD/$(basename $(basename ${outputFastq1} .fastq.gz) .fq.gz)/$(basename ${inputFastq2} )"_trimming_report.txt" ${outputFastq2}"_trimming_report.txt"
        fi
    }

    output {
        File fastq1 = outputFastq1
        File fastq1Log = outputFastq1 + "_trimming_report.txt"
        File? fastq2 = outputFastq2
        File? fastq2Log = outputFastq2 + "_trimming_report.txt"
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
    }
}