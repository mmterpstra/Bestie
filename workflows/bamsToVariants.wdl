version 1.0

import "../structs.wdl" as structs
import "../tasks/common.wdl" as common
import "../tasks/trimgalore.wdl" as trimgalore
import "../tasks/cutadapt.wdl" as cutadapt
import "../tasks/fastqc.wdl" as fastqc
import "../tasks/picard.wdl" as picard
import "../tasks/fgbio.wdl" as fgbio
import "../tasks/alignment.wdl" as align
import "../tasks/gatk.wdl" as gatk
import "../tasks/gatk_mutect2.wdl" as mutect
import "../tasks/bcftools.wdl" as bcftools


import "../tasks/ichorcna.wdl" as ichorcna
import "../workflows/qc.wdl" as qc

workflow BamsToVariants {
    input {
        String picardModule = "picard/2.26.10-Java-8-LTS"
        String gatkModule = "GATK/4.2.4.1-Java-8-LTS"
        String samtoolsModule = "SAMtools/1.15.1-GCC-11.3.0"
        String bcftoolsModule = "BCFtools/1.21-GCC-12.2.0"
        Reference reference
        IndexedFile dbsnp
        #IndexedFile cosmic
        IndexedFile gnomadOnlyAfVcf
        IndexedFile panelOfNormalsVcf
        Array[IndexedFile] knownSites
        #Array[SampleDescriptor] sample
        File? sampleJson
        SampleConfig? sampleConfigIn 
        File targetIntervalList
        Int targetScatter
    }
    #    SampleConfig sampleConfig = read_json(sampleJson)
    SampleConfig sampleConfig = if (defined(sampleConfigIn)) then select_first([sampleConfigIn]) else read_json(select_first([sampleJson]))

    call picard.SplitAndPadIntervals as splitIntervals {
        input:
            picardModule = picardModule,
            inputIntervalListFile = targetIntervalList,
            outputPrefix = "Targets_padded_and_split",
            targetScatter = targetScatter,
            padding = 200
    }
    
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
    }
    call common.UniqueArray as uniqueNormals {
        input:
            array = select_all(sampleNormals)
    }

    scatter (scatteredtargetsIdx in range(length(splitIntervals.paddedScatteredIntervalList))) {
        call mutect.MuTect2JointCalling as mutect2 {
            input:
                gatkModule = gatkModule,
                reference = reference,
                germlineAfVcf = gnomadOnlyAfVcf,
                panelOfNormalsVcf = panelOfNormalsVcf,
                inputBams = gatherBams.link,
                inputBamIndexes = gatherBais.link,
                normals = uniqueNormals.outArray,
                targetIntervalList = splitIntervals.paddedScatteredIntervalList[scatteredtargetsIdx],
                outputVcfBasename = "mutect_joint_calls_scat"+scatteredtargetsIdx
        }

    }
    call mutect.LearnReadOrientationModel as learnReadOrientationModel {
        input:
            gatkModule = gatkModule,
            f1r2TarGz = mutect2.f1r2TarGz,
            outputArtifactPriorTableBasename = "read-orientation-model"
    }
    call mutect.MergeMutectStats as mergeMutectStats {
        input:
            gatkModule = gatkModule,
            inputMutectStats = mutect2.stats,
            outputMergedStats = "merged.stats"
    }
    scatter (scatteredtargetsIdx in range(length(mutect2.vcfOut))) {

        call mutect.FilterMutect as FilterMutectCalls {
            input:
                gatkModule = gatkModule,
                reference = reference,

                inputSomaticVcf = mutect2.vcfOut[scatteredtargetsIdx],
                stats = mergeMutectStats.stats,
                artifactPriorsTarGz=learnReadOrientationModel.artifactpriortable,
                outputVcfBasename="mutect_filterd_scat_"+scatteredtargetsIdx,
                targetIntervalList = splitIntervals.paddedScatteredIntervalList[scatteredtargetsIdx]
        }
    }
    scatter (sample in sampleConfig.samples) {
        IndexedFile sampleIndexedBamHc = select_first([sample.alignedReads])

        scatter (scatteredtargetsIdx in range(length(splitIntervals.paddedScatteredIntervalList))) {
        
            #haplotypecallergvcf
            call gatk.HaplotypeCallerGVcf as haplotypeCallerGvcf {
                input:
                    gatkModule = gatkModule,
                    reference = reference,
                    inputBam = sampleIndexedBamHc,
                    outputVcfBasename = sample.name + ".idx_" + scatteredtargetsIdx,
                    targetIntervalList = splitIntervals.paddedScatteredIntervalList[scatteredtargetsIdx]
            }
            
            #mutect
             
            #if (defined(sample.control) ){
            #    scatter (sampleControl in sampleConfig.samples) {
            #        if(sample.control == sampleControl.name){
            #            call mutect.MuTect2JointCalling as mutect2 {
            #                input:
            #                    gatkModule = gatkModule,
            #                    reference = reference,
            #                    germlineAfVcf = gnomadOnlyAfVcf,
            #                    panelOfNormalsVcf = panelOfNormalsVcf,
            #                    inputBams = sampleIndexedBam.file,
            #                    inputBamIndexes = sampleIndexedBam.index,
            #                    inputIntervalListFile = splitIntervals.paddedScatteredIntervalList[scatteredtargetsIdx]
            #            }
            #        }
            #    }
            #}
            #Lofreq
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
    call picard.GatherVcfs as gatherMutect {
        input:
            picardModule = picardModule,
            inputVcfs=FilterMutectCalls.vcf,
            outputPrefix="project_mutect2",
    }
    call bcftools.Index as indexMutect {
        input:
            bcftoolsModule = bcftoolsModule,
            inputVcf = gatherMutect.outputVcf
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
        #output files of workflow
        #needs to be debugged probably
        
        #fastqtobam certain output
        Array[IndexedFile] mutect2Vars = FilterMutectCalls.vcfOut
        IndexedFile mutect2Vcf = indexMutect.vcfOut
        IndexedFile haplotypecallergVcf = gatherHcSamples.vcfOut
        IndexedFile haplotypecallerVcf = genotypeHcProjectGvcf.vcfOut
        
    }

    meta {
        author: "MMTerpstra"
        description: "This is the single samples bams to variants workflow."
    }  
}