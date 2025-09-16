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
import "../tasks/lofreq.wdl" as lofreq
import "../tasks/pipeline-util.wdl" as util
import "../tasks/ichorcna.wdl" as ichorcna
import "../workflows/qc.wdl" as qc
import "../workflows/bamsToGermlineVariants.wdl" as bamsToGermlineVariants


workflow BamsToSomaticMutectVariants {
    input {
        String picardModule = "picard/2.26.10-Java-8-LTS"
        String gatkModule = "GATK/4.2.4.1-Java-8-LTS"
        String bcftoolsModule = "BCFtools/1.21-GCC-12.2.0"
        String pipelineUtilModule = "pipeline-util/0.8.20-5-ga0a29bb-foss-2024a"
        Reference reference
        IndexedFile dbsnp
        IndexedFile gnomadOnlyAfVcf
        IndexedFile panelOfNormalsVcf
        Array[IndexedFile] knownSites
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
    #Array[String] sampleNames 
    call common.UniqueArray as uniqueNormals {
        input:
            array = select_all(sampleNormals)
    }
    #this should be the best way of selecting normal samples 
    #idk how it goes for wdl validation though
    #Array[SampleDescriptor] normalSamples = []

    #this works under the assumtion that at least one sampledescriptor is present for each normal name
    
    #better to post filter on normal samples 
    #scatter (sample in viterbiSample) {
    #    scatter(normalname in select_all(uniqueNormals.outArray)){
    #        if (sample.name == normalname){
    #            Boolean isTumor = false 
    #        }
    #    }
    #    Boolean sampleIsTumor = select_first(flatten([isTumor,[true]]))
    #    placeholder for dbugging
    #    call common.UniqueArray as uniqueInputs {
    #        input:
    #            array = select_all(flatten([[sampleIsTumor],[isTumor,[true]]]))
    #    }
    #    if(sampleIsTumor){
    #        SampleDescriptor viterbiTumorSampleScat = sample
    #    }
    #}
    #Array[SampleDescriptor] viterbiTumor = select_all(viterbiTumorSampleScat)


    #this might also be a way of solving the upper problem of merging the optional normalMatched samples to a list of normal samples
    #call common.SelectNormals as selectNormalSamples {
    #    input:
    #        normaleNames = select_all(uniqueNormals.outArray),
    #        samples = sampleConfig.samples
    #}

    scatter (scatteredtargetsIdx in range(length(targetIntervals))) {
        call mutect.MuTect2JointCalling as mutect2 {
            input:
                gatkModule = gatkModule,
                reference = reference,
                germlineAfVcf = gnomadOnlyAfVcf,
                panelOfNormalsVcf = panelOfNormalsVcf,
                inputBams = gatherBams.link,
                inputBamIndexes = gatherBais.link,
                normals = uniqueNormals.outArray,
                targetIntervalList = targetIntervals[scatteredtargetsIdx],
                outputVcfBasename = "mutect_joint_calls_scat"+scatteredtargetsIdx,
                targetScatter = length(targetIntervals)
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
    scatter (scatteredtargetsIdx in range(length(targetIntervals))) {
        #Mutect downstream
        call mutect.FilterMutect as filterMutectCalls {
            input:
                gatkModule = gatkModule,
                reference = reference,
                inputSomaticVcf = mutect2.vcfOut[scatteredtargetsIdx],
                stats = mergeMutectStats.stats,
                artifactPriorsTarGz=learnReadOrientationModel.artifactpriortable,
                outputVcfBasename="mutect_filterd_scat_"+scatteredtargetsIdx,
                targetIntervalList = targetIntervals[scatteredtargetsIdx]
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
    }
    
    call bcftools.Stats as mutectRawStats {
        input:
            inputVcfs=mutect2.vcf,
            inputIndexedVcfs=mutect2.vcfOut,
            outputBasename="mutect_raw",
            bcftoolsModule=bcftoolsModule,
    }

    call bcftools.Stats as mutectFilteredStats {
        input:
            inputVcfs=normaliseMutect.vcf,
            inputIndexedVcfs=normaliseMutect.vcfOut,
            outputBasename="mutect_filtered",
            bcftoolsModule=bcftoolsModule,
    }

    #merge all regions

    output {
        #output files of workflow
        #needs to be debugged probably
        Array[File] vcf = normaliseMutect.vcf
        Array[IndexedFile] idxVcf = normaliseMutect.vcfOut
        Array[File] bcftoolsStats = [mutectRawStats.stats,mutectFilteredStats.stats]
    }

    meta {
        author: "MMTerpstra"
        description: "This is the single samples bams to variants workflow."
    }  
}