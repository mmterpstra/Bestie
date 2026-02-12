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
        Int timeMinutes = 1 + ceil(size(inputFastq1, "G")) * 30
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

task CutadaptUbam {
    input {
        File inputUbam
        String outputFastq1
        String? outputFastq2
        Int minimumLength = 15
        Int baseQuality = 10
        Int minAdapterOverlap = 1
        Int threads = 1
        Array[String] read1Adapters = ["AGATCGGAAGAGC"]
        Array[String] read2Adapters = ["AGATCGGAAGAGC"]
        Int? memoryGb = 4
        Int javaXmxMemoryMb = floor((memoryGb-2)*0.95*1024)
        Int timeMinutes = 1 + ceil(size(inputUbam, "G")) * 50
        Int disk = ceil(size(inputUbam, "M")*2.1)

        #File? fastq_input_umi
        #String samplename
        #String identifier
	    String cutadaptModule
        String picardModule
    }
    command {
        set -e 
        module load ~{picardModule}

        echo "time: "~{timeMinutes}";disk:"~{disk}
        #intermediate output file removal or storage on $TMP_DIR
        java -Xmx~{javaXmxMemoryMb}m \
            -jar $EBROOTPICARD/picard.jar SamToFastq \
            INPUT=~{inputUbam} \
            OUTPUT_PER_RG=false \
            --FASTQ ~{outputFastq1}.tmp_1.fastq.gz \
            ~{ "--SECOND_END_FASTQ " + outputFastq2 + ".tmp_2.fastq.gz"}

        module load ${cutadaptModule}

        if [ -e ~{outputFastq2}".tmp_2.fastq.gz" ];then \
            cutadapt \
                --minimum-length ~{minimumLength} \
                -j ~{threads} \
                -e 0.1 \
                -q ~{baseQuality} \
                -O ~{minAdapterOverlap} \
                -a ~{sep=' -a ' read1Adapters} \
                --output ~{outputFastq1} \
                ~{outputFastq1}.tmp_1.fastq.gz \
                &>  ~{outputFastq1}"_trimming_report.txt"
        else \
            cutadapt \
                --minimum-length ~{minimumLength} \
                -j ~{threads} \
                -e 0.1 \
                -q ~{baseQuality} \
                -O ~{minAdapterOverlap} \
                -a ~{sep=' -a ' read1Adapters} \
                -A ~{sep=' -A ' read2Adapters} \
                --output ~{outputFastq1} \
                --paired-output ~{outputFastq2} \
                ~{outputFastq1}.tmp_1.fastq.gz \
                ~{outputFastq2}.tmp_2.fastq.gz \
                &>  ~{outputFastq1}"_trimming_report.txt"

                rm -v ~{outputFastq2}.tmp_2.fastq.gz
        fi

        rm -v ~{outputFastq1}.tmp_1.fastq.gz
    }

    output {
        File fastq1 = outputFastq1
        File fastq1Log = outputFastq1 + "_trimming_report.txt"
        File? fastq2 = outputFastq2
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
    }
}