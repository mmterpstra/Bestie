version 1.0

import "../structs.wdl" as structs
import "../tasks/common.wdl" as common
import "../tasks/gatk.wdl" as gatk
import "../tasks/gatk_mutect2.wdl" as mutect

workflow BamsToGermlineVariants {
    input {
        String gatkModule = "GATK/4.2.4.1-Java-8-LTS"
        Reference reference
        File? sampleJson
        SampleConfig? sampleConfigIn 
        Array[File] targetIntervals
    }
    #This is a toggle to either use the sampleconfig in the input json or as separate file for conveniance 
    SampleConfig sampleConfig = if (defined(sampleConfigIn)) then select_first([sampleConfigIn]) else read_json(select_first([sampleJson]))

    
    #Array[IndexedFile] sampleIndexedBams = sampleConfig[].samples.alignedReads
    scatter (sample in sampleConfig.samples) {
        IndexedFile sampleIndexedBam = select_first([sample.alignedReads])

        call common.CreateLink as gatherBams {
            input:
                inputFile = sampleIndexedBam.file,
                hardLink = true,
                outputPath = select_first([sample.name])+ "_aligned_reads.bam",
        }
        call common.CreateLink as gatherBais {
            input:
                inputFile = sampleIndexedBam.index,
                hardLink = true,
                outputPath = select_first([sample.name])+ "_aligned_reads.bai",
        }
        if(defined(sample.normal)){
            String sampleNormals = select_first([sample.normal])
        }
        String sampleNames = sample.name
    }

    #haplotypecaller gvcf mode workflow
    scatter (sample in sampleConfig.samples) {
        IndexedFile sampleIndexedBamHc = select_first([sample.alignedReads])

        scatter (scatteredtargetsIdx in range(length(targetIntervals))) {
        
            #haplotypecallergvcf
            call gatk.HaplotypeCallerGVcf as haplotypeCallerGvcf {
                input:
                    gatkModule = gatkModule,
                    reference = reference,
                    inputBam = sampleIndexedBamHc,
                    outputVcfBasename = sample.name + ".idx_" + scatteredtargetsIdx,
                    targetIntervalList = targetIntervals[scatteredtargetsIdx],
                    targetScatter=length(targetIntervals)
            }
        }
        #CombineGVCFs to create a single sample gvcf
        call gatk.CombineGVCFs as gatherRegions {
                input:
                    gatkModule = gatkModule,
                    reference = reference,
                    inputGVcfs = haplotypeCallerGvcf.vcfOut,
                    inputGVcfsFiles = haplotypeCallerGvcf.vcf,
                    outputVcfBasename = sample.name + "_hcgvcf",
        }

    }

    #CombineGVCFs to create a single project gvcf
    call gatk.CombineGVCFs as gatherHcSamples {
        input:
            gatkModule = gatkModule,
            reference = reference,
            inputGVcfs = gatherRegions.vcfOut,
            inputGVcfsFiles = gatherRegions.vcf,
            outputVcfBasename = "hcgvcf_allsamples",
    }
    #genotype gvcf
    call gatk.GenotypeGVCFs as genotypeHcProjectGvcf {
        input:
            gatkModule = gatkModule,
            reference = reference,
            inputGVcf = gatherHcSamples.vcfOut,
            inputGVcfsFile = gatherHcSamples.vcf,
            outputVcfBasename = "hc_allsamples",
    }

    output {
        IndexedFile idxGvcf = gatherHcSamples.vcfOut
        IndexedFile idxVcf = genotypeHcProjectGvcf.vcfOut
        #Array[IndexedFile] = 
        #Array[File] = 
    }

    meta {
        author: "MMTerpstra"
        description: "This is the single samples bams to variants workflow."
    }  
}