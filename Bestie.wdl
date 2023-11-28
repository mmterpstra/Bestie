version 1.0

import "structs.wdl" as structs
import "tasks/common.wdl" as common
import "tasks/trimgalore.wdl" as trimgalore
import "tasks/cutadapt.wdl" as cutadapt
import "tasks/fastqc.wdl" as fastqc
import "tasks/picard.wdl" as picard
import "tasks/fgbio.wdl" as fgbio
import "tasks/alignment.wdl" as align
import "tasks/gatk.wdl" as gatk
import "tasks/ichorcna.wdl" as ichorcna

import "workflows/qc.wdl" as qc
import "workflows/fastqToBam.wdl" as fastqToBam

import "tasks/multiqc.wdl" as multiqc


workflow FastqToVariants {
    input {
        String fastqcModule = "FastQC/0.11.9-Java-11"
        String trimgaloreModule = "Trim_Galore/0.6.6-GCCcore-9.3.0-Python-3.8.2"
        String multiqcModule = "multiqc/1.12-GCCcore-11.3.0"
        String bwaModule = "BWA/0.7.17-GCCcore-11.3.0"
        String picardModule = "picard/2.26.10-Java-8-LTS"
        String gatkModule = "GATK/4.2.4.1-Java-8-LTS"
        String hmmcopyutilsModule = "hmmcopy_utils/5911bf69f1-foss-2022a"
        String samtoolsModule = "SAMtools/1.15.1-GCC-11.3.0"
        String fgbioModule = "fgbio/1.3.0"
        Boolean runCutadapt = false 
        String cutadaptModule = "cutadapt/4.2-GCCcore-11.3.0"
        Array[String] read1Adapters = ["AGATCGGAAGAGC"]
        Array[String] read2Adapters = ["AGATCGGAAGAGC"]
        Boolean runTwistUmi = false
        Boolean runReadcounter = true
        Boolean runBaseQualityRecalibration = true
        File sampleJson
        Reference reference
        IndexedFile dbsnp
        IndexedFile cosmic
        Array[IndexedFile] knownSites
        File targetIntervalList
        Int targetScatter
        BwaIndex referenceBwaIndex
    }
    SampleConfig sampleConfig = read_json(sampleJson)


    call picard.SplitAndPadIntervals as splitIntervals {
    input:
        picardModule = picardModule,
        inputIntervalListFile = targetIntervalList,
        outputPrefix = "Targets_padded_and_split",
        targetScatter = targetScatter,
        padding = 200
    }
    scatter (sample in sampleConfig.samples) {
        call fastqToBam.FastqToBam as fqToBam {
            input:
                sample = sample,
                fastqcModule = fastqcModule,
                trimgaloreModule = trimgaloreModule,
                bwaModule = bwaModule,
                picardModule = picardModule,
                gatkModule = gatkModule,
                samtoolsModule = samtoolsModule,
                fgbioModule = fgbioModule,
                runCutadapt = runCutadapt,
                cutadaptModule = cutadaptModule,
                read1Adapters = read1Adapters,
                read2Adapters = read2Adapters,
                runBaseQualityRecalibration = runBaseQualityRecalibration,
                reference = reference,
                referenceBwaIndex = referenceBwaIndex,
                runCutadapt = runCutadapt,
                dbsnp = dbsnp,
                knownSites = knownSites,
                targetIntervalList = targetIntervalList
        }
        #Does not work: sample.alignedReads = fqToBam. That is why I added this task to serialise/deserialise the json object and add the variables.
        call common.AddAlignedReadsToSampleDescriptor as addBamToSample {
            input:
                sample=sample,
                bam=fqToBam.bam
        }
        
        SampleDescriptor sampleNew = addBamToSample.sampleUpdated

        if(runReadcounter){
            call ichorcna.HmmcopyReadcounter as readcounter500kbp {
                input:
                    inputBam=fqToBam.bam.file,
                    inputBai=fqToBam.bam.index,
                    outputPrefix=sample.name,
                    windowkilobase=500,
                    referencefai=reference.fai,
                    hmmcopyutilsModule=hmmcopyutilsModule,
                    samtoolsModule=samtoolsModule
            }
            call ichorcna.HmmcopyReadcounter as readcounter1000kbp {
                input:
                    inputBam=fqToBam.bam.file,
                    inputBai=fqToBam.bam.index,
                    outputPrefix=sample.name,
                    windowkilobase=1000,
                    referencefai=reference.fai,
                    hmmcopyutilsModule=hmmcopyutilsModule,
                    samtoolsModule=samtoolsModule
            }
        }
        
        #variant callin to sub pipeline

        #scatter by sequencing targets intervals
        #scatter (scatteredtargetsIdx in range(length(splitIntervals.paddedScatteredIntervalList))) {
            #haplotypecallergvcf
        #    call gatk.HaplotypeCallerGVcf as haplotypeCallerGvcf {
        #        input:
        #           gatkModule = gatkModule,
        #            reference = reference,
        #            inputBam = fqToBam.bam,
        #            outputVcfBasename = sample.name + ".idx_" + scatteredtargetsIdx,
        #            targetIntervalList = splitIntervals.paddedScatteredIntervalList[scatteredtargetsIdx]
        #    }
            
            #genotypegvcf per sample?

            #mutect

            #Lofreq
        #}

    }
    #run multiqc to bundle outputs
    Array[File] files = flatten(
        flatten(
            [
                fqToBam.fastqcZip,
                [
                    fqToBam.markdupLog,
                    fqToBam.qcZip,
                ]
            ]
        )
    )
    #Array[String] optionalFiles = 
    call multiqc.MultiQC as multiqc {
            input:
                multiqcModule = multiqcModule,
                files = files,
                optionalFiles = select_all(flatten(flatten([fqToBam.cutadaptLogs,[fqToBam.umiQcZip]])))
    } 
    output {
        #output files of workflow
        #needs to be debugged probably
        
        #fastqtobam certain output
        Array [File] markdupLogs = fqToBam.markdupLog
        Array [File] qcZips = fqToBam.qcZip
        #fastqtobam optional output
        Array [File?] cutadaptLogs = flatten(fqToBam.cutadaptLogs)
        Array [File?] umiQcZips = fqToBam.umiQcZip
        Array [IndexedFile] bams = fqToBam.bam

        Array[SampleDescriptor] samples = SampleNew
        #This depends on fastqToBam
        File multiqcHtml = multiqc.html
        File multiqcDir = multiqc.dir
    }

    meta {
        author: "MMTerpstra"
        description: "This is the single samples fastqs to aligned bam workflow."
    }  
}