version 1.0

task MarkTrimming {
    input {
        #query sorted unaligned bam
        File inputUbam
        File cutadaptFastq1
        File? cutadaptFastq2
        String outputTrimmedBamBase
        Int memoryGb = "1"
        #Int javaXmxMemoryMb = floor((memoryGb-0.5)*0.95*1024)
        String marktrimmingModule = "marktrimming/0.0.2-GCC-12.2.0"
        Int timeMinutes = 5 + ceil(size(select_all([inputUbam,cutadaptFastq1,cutadaptFastq2]), "G")) * 50
    }
    command {
        #ideally in this kind of a workflow
        #(ml Pysam ;cat unaligned.bam| python3 ../marktrimming/build/scripts-3.10/marktrimming.py --ubam-out True --output - --input - --fastq cutadapt_R1.fastq.gz --fastq cutadapt_R2.fastq.gz )| (ml SAMtools ;samtools view -h - |perl -wpe 's/XT:i:123/XT:i:1/')|(ml SAMtools ; samtools view -Sb - )| (ml picard && java -jar $EBROOTPICARD/picard.jar SamToFastq I=/dev/stdin INTERLEAVE=true CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 FASTQ=/dev/stdout )
        set -e
        module load ~{marktrimmingModule} && \
        marktrimming \
            --input ~{inputUbam} \
            --output ~{outputTrimmedBamBase}'.bam' \
            --fastq ~{cutadaptFastq1} \
            ~{"--fastq " + cutadaptFastq2}   
    }

    output {
        File bam = outputTrimmedBamBase + '.bam'
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
    }
}