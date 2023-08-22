version 1.0

import "structs.wdl" as structs
import "tasks/common.wdl" as common
import "tasks/trimgalore.wdl" as trimgalore
import "tasks/fastqc.wdl" as fastqc

import "tasks/picard.wdl" as picard
import "tasks/alignment.wdl" as align
import "tasks/gatk.wdl" as gatk

#import "tasks/ichorcna.wdl" as ichorcna

import "workflows/qc.wdl" as qc


import "tasks/multiqc.wdl" as multiqc


workflow fastqToVariants {
    input {
        String fastqcModule = "FastQC/0.11.9-Java-11"
        String trimgaloreModule = "Trim_Galore/0.6.6-GCCcore-9.3.0-Python-3.8.2"
        String multiqcModule = "multiqc/1.12-GCCcore-11.3.0"
        String bwaModule = "BWA/0.7.17-GCCcore-11.3.0"
        String picardModule = "picard/2.26.10-Java-8-LTS"
        String gatkModule = "GATK/4.2.4.1-Java-8-LTS"
        String multiqcModule = "multiqc/1.12-GCCcore-11.3.0"
        String hmmcopyutilsModule = "hmmcopy_utils/5911bf69f1-foss-2022a"
        String samtoolsModule = "SAMtools/1.15.1-GCC-11.3.0"
        Boolean runReadcounter = true
        File sampleJson
        Reference reference
        File targetIntervalList
        Int targetScatter
        BwaIndex referenceBwaIndex
        #Array[File] = knownIndelsVcfs
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

        scatter (rg in sample.readgroups) {
            #linking for uniform filenames
            call common.CreateLink as getfastq1 {
            input:
                inputFile = rg.fastq1,
                outputPath = sample.name + "_" + rg.flowcell + "_" + rg.identifier + "_R1.fastq.gz"
            }
            call fastqc.FastQC as fastqc1 {
            input:
                inputFastq = getfastq1.link,
                fastqcModule = fastqcModule
            }
            if (defined(rg.fastq2)) {
                call common.CreateLink as getfastq2 {
                input:
                    inputFile = select_first([rg.fastq2]),
                    outputPath = sample.name + "_" + rg.flowcell + "_" + rg.identifier + "_R2.fastq.gz"
                }
                call fastqc.FastQC as fastqc2 {
                input:
                    inputFastq = getfastq2.link,
                    fastqcModule = fastqcModule
                }

            }
            #trim adapters
            call trimgalore.TrimGalore as adaptertrim {
                input:
                    inputFastq1 = getfastq1.link,
                    inputFastq2 = getfastq2.link,
                    outputFastq1 = sample.name + "_" + rg.identifier + "_trim_R1.fastq.gz",
                    outputFastq2 = sample.name + "_" + rg.identifier + "_trim_R2.fastq.gz",
                    memoryGb = 1,
                    trimgaloreModule = trimgaloreModule
            }
            
            #align
            ##to samconversion
            call picard.FastqToUnmappedBam as fastqToBam {
                input:
                    inputFastq1 = adaptertrim.fastq1,
                    inputFastq2 = adaptertrim.fastq2,
                    picardModule = picardModule,
                    sampleName = sample.name, 
                    platform = rg.platform,
                    platformUnit = rg.run  + "_" + rg.barcode1 + "+" + select_first([rg.barcode2,'AAAAAA']) + "." + rg.lane,
                    readGroupName = rg.run  + "_" + rg.flowcell  + "_" + rg.barcode1 + "+" + select_first([rg.barcode2,'AAAAAA']) + "." + rg.lane,
                    outputUnalignedBam = rg.run  + "_" + rg.flowcell  + "_" + rg.barcode1 + "+" + select_first([rg.barcode2,'AAAAAA']) + "." + rg.lane + "_unaligned.bam",
            }
            ##map with bwa
            call align.bwaAlignBam as bwaBam {
                input:
                    inputUnalignedBam = fastqToBam.unalignedBam,
                    referenceBwaIndex = referenceBwaIndex,
                    reference = reference,
                    bwaModule = bwaModule,
                    picardModule = picardModule,
                    outputBam = rg.run  + "_" + rg.flowcell  + "_" + rg.barcode1 + "+" + select_first([rg.barcode2,'AAAAAA']) + "." + rg.lane + "_aligned.bam",
            }
        }

        call fastqc.FastQCSample as fastqcSample1 {
            input:
                fastqcModule = fastqcModule,
                inputFastqGzs = getfastq1.link,
                outputPrefix = sample.name + "_R1",
        }
        if (defined(getfastq2.link)) {
            call fastqc.FastQCSample as fastqcSample2 {
            input:
                fastqcModule = fastqcModule,
                inputFastqGzs = select_all(getfastq2.link),
                outputPrefix = sample.name + "_R2",
            }
        }

        call picard.MarkDuplicates as markDups {
            input:
                picardModule = picardModule,
                inputBams = bwaBam.alignedBam,
                outputBamBasename = sample.name + '_markdup',
                outputMetrics = sample.name + '.markdup_metrics'

        }
        call picard.SortSam as sortBam {
            input: 
                picardModule = picardModule,
                inputBam = markDups.bam,
                outputBamBasename = sample.name + '_sort'
                
        }
        #if(runReadcounter){
        #    call ichorcna.hmmcopyReadcounter as readcounter {
        #        input:
        #            inputBam=sortBam.bam,
        #            outputPrefix=sample.name,
        #            windowkilobase=500,
        #            hmmcopyutilsModule=hmmcopyutilsModule,
        #            samtoolsModule=samtoolsModule
        #    }
        #}
        #optional bqsr
        #optional indelrealignment

        #run qc
        call qc.bamQualityControl as bamQualityControl {
        #call gatk.CollectMultipleMetrics as CollectMultipleMetrics {
            input:
                gatkModule = gatkModule,
                picardModule = picardModule,
                reference = reference,
                inputBam = sortBam.bam,
                inputBai = sortBam.bai,
                outputPrefix =  sample.name + '_qc',
                targetIntervalList = targetIntervalList
        }
        #scatter by sequencing targets intervals
        scatter (scatteredtargetsIdx in range(length(splitIntervals.paddedScatteredIntervalList))) {
            #haplotypecallergvcf
            call gatk.HaplotypeCaller as haplotypeCallerGvcf {
                input:
                    gatkModule = gatkModule,
                    reference = reference,
                    inputBam = sortBam.bam,
                    inputBai = sortBam.bai,
                    outputVcfBasename = sample.name + ".idx_" + scatteredtargetsIdx,
                    targetIntervalList = splitIntervals.paddedScatteredIntervalList[scatteredtargetsIdx]
            }
            
            #genotypegvcf

            #mutect

            #Lofreq
        }
    }

    #run multiqc to bundle outputs
    Array[File] files = flatten(
        flatten(
            [
                fastqc1.outZip,adaptertrim.fastq1Log,
                [
                    markDups.metrics,
                    bamQualityControl.alignmentMetrics,
                    bamQualityControl.baseDistributionMetrics,
                    bamQualityControl.baseDistributionPdf,
                    bamQualityControl.insertSizeMetrics,
                    bamQualityControl.insertSizePdf,
                    bamQualityControl.qualityByCycleMetrics,
                    bamQualityControl.qualityByCyclePdf,
                    bamQualityControl.qualityDistributionMetrics,
                    bamQualityControl.qualityDistributionPdf,
                    bamQualityControl.readLengthPdf,
                    fastqcSample1.outZip
                ]
            ]
        )
    )
    #Array[String] optionalFiles = 
    call multiqc.MultiQC as multiqc {
            input:
                multiqcModule = multiqcModule,
                files = files,
                optionalFiles = select_all(flatten(flatten([fastqc2.outZip,[fastqcSample2.outZip, bamQualityControl.hsMetrics]])))
    }
    
    
}