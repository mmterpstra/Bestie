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
import "../tasks/freebayes.wdl" as freebayes
import "../tasks/pipeline-util.wdl" as util

import "../tasks/ichorcna.wdl" as ichorcna
import "../workflows/qc.wdl" as qc

workflow BamsToVariants {
    input {
        String picardModule = "picard/2.26.10-Java-8-LTS"
        String gatkModule = "GATK/4.2.4.1-Java-8-LTS"
        String samtoolsModule = "SAMtools/1.15.1-GCC-11.3.0"
        String bcftoolsModule = "BCFtools/1.21-GCC-12.2.0"
        String freebayesModule = "freebayes/1.3.7-gfbf-2024a-R-4.4.2"
        String pipelineUtilModule = "pipeline-util/0.8.20-5-ga0a29bb-foss-2024a"
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
        call freebayes.FreebayesSomatic as freebayesSomatic {
            input:
                freebayesModule = freebayesModule,
                reference = reference,
                inputBams = gatherBams.link,
                inputBamIndexes = gatherBais.link,
                targetIntervalList = splitIntervals.paddedScatteredIntervalList[scatteredtargetsIdx],
                outputVcfBasename = "freebayes_joint_calls_scat"+scatteredtargetsIdx
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
    scatter (scatteredtargetsIdx in range(length(splitIntervals.paddedScatteredIntervalList))) {
        #Mutect downstream
        call mutect.FilterMutect as filterMutectCalls {
            input:
                gatkModule = gatkModule,
                reference = reference,
                inputSomaticVcf = mutect2.vcfOut[scatteredtargetsIdx],
                stats = mergeMutectStats.stats,
                artifactPriorsTarGz=learnReadOrientationModel.artifactpriortable,
                outputVcfBasename="mutect_filterd_scat_"+scatteredtargetsIdx,
                targetIntervalList = splitIntervals.paddedScatteredIntervalList[scatteredtargetsIdx]
        }
        call util.Callerise as mutectScatCallerise {
            input:
                pipelineUtilModule=pipelineUtilModule,
                inputVcf = filterMutectCalls.vcf,
                outputBase = "mutect_tagged_scat_"+scatteredtargetsIdx,
                caller="MuTect2_"
        }
        call bcftools.Norm as normaliseMutect {
            input:
                bcftoolsModule=bcftoolsModule,
                reference=reference,
                inputVcf=mutectScatCallerise.vcf,
                outputBasename ="mutect_norm_scat_"+scatteredtargetsIdx,
        }
        #freebayes downstream
        call bcftools.ViewSamples as freebayesNormals {
            input:
                bcftoolsModule=bcftoolsModule,
                inputVcf = freebayesSomatic.vcf[scatteredtargetsIdx],
                outputBasename = "freebayes_normals_scat_"+scatteredtargetsIdx,
                samples=uniqueNormals.outArray
        }
        call util.AdFilter as freebayesFilterGermline {
            input:
                inputVariantsToFilter = freebayesSomatic.vcfOut[scatteredtargetsIdx],
                inputVcfsFiles=[freebayesNormals.vcf],
                inputVcfs=[freebayesNormals.vcfOut],
                outputBasename='freebayes_normals_filt_scat_'+ scatteredtargetsIdx,
        }
        #filterfreebayes?
        call util.Callerise as freebayesScatCallerise {
            input:
                pipelineUtilModule=pipelineUtilModule,
                inputVcf = freebayesFilterGermline.vcf,
                outputBase = "freebayes_tagged_scat_"+scatteredtargetsIdx,
                caller="freebayes_"
        }
        call bcftools.Norm as normaliseFreeBayes {
            input:
                bcftoolsModule=bcftoolsModule,
                reference=reference,
                inputVcf=freebayesScatCallerise.vcf,
                outputBasename="freebayes_norm_scat_"+scatteredtargetsIdx,
        }
        #call gatk.GenomicsDBImport as importCallers {
        #    input:
        #        gatkModule = gatkModule,
        #        reference = reference,
        #        inputVcfsFiles=[normaliseFreeBayes.vcf,normaliseMutect.vcf],
        #        inputVcfs=[normaliseFreeBayes.vcfOut,normaliseMutect.vcfOut],
        #        outputBasename='variants_merged_scat'+ +scatteredtargetsIdx,
        #        interval_list=splitIntervals.paddedScatteredIntervalList[scatteredtargetsIdx],
        #}
        #call gatk.SelectVariants as mergedVariantVcf {
        #    input:
        #        gatkModule = gatkModule,
        #        reference = reference,
        #        genomicsDbTar = importCallers.genomicsDbTar,
        #        outputVcfBasename = 'merged_scat'+ +scatteredtargetsIdx,
        #        
        #
        #}
        call bcftools.Isec as mergeVariantLists {
            input:
                bcftoolsModule=bcftoolsModule,
                inputVcfsFiles=[normaliseFreeBayes.vcf,normaliseMutect.vcf],
                inputVcfs=[normaliseFreeBayes.vcfOut,normaliseMutect.vcfOut],
                outputBasename='variants_merged_scat'+ scatteredtargetsIdx,
        }
        call freebayes.FreebayesRecall as annotateMergedLists {
            input:
                freebayesModule = freebayesModule,
                reference = reference,
                inputBams = gatherBams.link,
                inputBamIndexes = gatherBais.link,
                inputVariants=mergeVariantLists.vcfOut,
                outputVcfBasename = "variants_freebayesanno_merged_scat"+scatteredtargetsIdx
        }
        call bcftools.Norm as normaliseMerged {
            input:
                bcftoolsModule=bcftoolsModule,
                reference=reference,
                inputVcf=annotateMergedLists.vcfOut.file,
                outputBasename ="norm_merged_scat_"+scatteredtargetsIdx,
        }
        call util.ReannotateVariants as annotateCallers {
            input:
                pipelineUtilModule=pipelineUtilModule,
                combinedVariants=normaliseMerged.vcfOut,
                inputVcfsFiles=[normaliseFreeBayes.vcf,normaliseMutect.vcf],
                inputVcfs=[normaliseFreeBayes.vcfOut,normaliseMutect.vcfOut],
                outputBasename='annotCallers_merged_scat'+ scatteredtargetsIdx,
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
    call picard.SortVcfsIndexed as gatherMutect {
        input:
            picardModule = picardModule,
            inputVcfs=filterMutectCalls.vcf,
            outputPrefix="project_mutect2",
            createIndex=true
    }
    call picard.SortVcfsIndexed as gatherFreebayes {
        input:
            picardModule = picardModule,
            inputVcfs=freebayesSomatic.vcf,
            outputPrefix="project_mutect2",
            createIndex=true
    }
    #call bcftools.Index as indexMutect {
    #    input:
    #        bcftoolsModule = bcftoolsModule,
    #        inputVcf = gatherMutect.outputVcf
    #} 
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
        Array[IndexedFile] mutect2Vars = filterMutectCalls.vcfOut
        #IndexedFile mutect2Vcf = indexMutect.vcfOut
        IndexedFile haplotypecallergVcf = gatherHcSamples.vcfOut
        IndexedFile haplotypecallerVcf = genotypeHcProjectGvcf.vcfOut
        
    }

    meta {
        author: "MMTerpstra"
        description: "This is the single samples bams to variants workflow."
    }  
}