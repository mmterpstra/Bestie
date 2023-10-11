version 1.0

task Cutadapt {
    input {
        File inputFastq1
        File? inputFastq2
        String outputFastq1
        String? outputFastq2
        Int minimumLength = 15
        Array[String] read1Adapters = ["AGATCGGAAGAGC"]
        Array[String] read2Adapters = ["AGATCGGAAGAGC"]
        Int? memoryGb = 1
        Int timeMinutes = 1 + ceil(size(inputFastq1, "G")) * 50
        #File? fastq_input_umi
        #String samplename
        #String identifier
	    String cutadaptModule
    }
    command {
        set -e 
        module load ${cutadaptModule} && \
        if [ "${inputFastq2}x" == "x" ];then \
            cutadapt \
                --minimum-length ~{minimumLength} \
                -j 1 -e 0.1 -q 20 -O 1 \
                -a ~{sep=' -a ' read1Adapters} \
                --output ~{outputFastq1} \
                ~{inputFastq1} \
                &>  ~{outputFastq1}"_trimming_report.txt"
            #cutadapt  -j 1 -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC R2.gz
            #trim_galore "${inputFastq1}" \
            #    --output_dir "$(basename $(basename ${outputFastq1} .fastq.gz) .fq.gz)"
            #    ln -sf $PWD/$(basename $(basename ${outputFastq1} .fastq.gz) .fq.gz)/$(basename $(basename ${inputFastq1} .fastq.gz) .fq.gz)_trimmed.fq.gz ${outputFastq1}
            #    ln -sf $PWD/$(basename $(basename ${outputFastq1} .fastq.gz) .fq.gz)/$(basename ${inputFastq1} )"_trimming_report.txt" ${outputFastq1}"_trimming_report.txt"
        else \
            cutadapt \
                --minimum-length ~{minimumLength} \
                -j 1 -e 0.1 -q 20 -O 1 \
                -a ~{sep=' -a ' read1Adapters} \
                -A ~{sep=' -A ' read2Adapters} \
                --output ~{outputFastq1} \
                --paired-output ~{outputFastq2} \
                ~{inputFastq1} \
                ~{inputFastq2} \
                &>  ~{outputFastq1}"_trimming_report.txt"
            #-a~{sep=' -a' read1Adapters}
            #trim_galore --paired "${inputFastq1}" "${inputFastq2}" \
            #  --output_dir "$(basename $(basename ${outputFastq1} .fastq.gz) .fq.gz)" 
            #ln -sf $PWD/$(basename $(basename ${outputFastq1} .fastq.gz) .fq.gz)/$(basename $(basename ${inputFastq1} .fastq.gz) .fq.gz)_val_1.fq.gz ${outputFastq1}
            #ln -sf $PWD/$(basename $(basename ${outputFastq1} .fastq.gz) .fq.gz)/$(basename ${inputFastq1} )"_trimming_report.txt" ${outputFastq1}"_trimming_report.txt"
            #ln -sf $PWD/$(basename $(basename ${outputFastq1} .fastq.gz) .fq.gz)/$(basename $(basename ${inputFastq2} .fastq.gz) .fq.gz)_val_2.fq.gz ${outputFastq2}
            #ln -sf $PWD/$(basename $(basename ${outputFastq1} .fastq.gz) .fq.gz)/$(basename ${inputFastq2} )"_trimming_report.txt" ${outputFastq2}"_trimming_report.txt"
        fi
    }

    output {
        File fastq1 = outputFastq1
        File fastq1Log = fastq1 + "_trimming_report.txt"
        File? fastq2 = outputFastq2
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
    }
}