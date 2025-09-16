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
import "../workflows/bamsToSomaticFreebayesVariants.wdl" as somaticFreebayes
import "../workflows/bamsToSomaticMutectVariants.wdl" as somaticMutect
import "../workflows/bamsToSomaticLofreqVariants.wdl" as somaticLofreq


workflow BamsToSomaticVariants {
    input {
        String picardModule = "picard/2.26.10-Java-8-LTS"
        String gatkModule = "GATK/4.2.4.1-Java-8-LTS"
        String samtoolsModule = "SAMtools/1.15.1-GCC-11.3.0"
        String bcftoolsModule = "BCFtools/1.21-GCC-12.2.0"
        String freebayesModule = "freebayes/1.3.7-gfbf-2024a-R-4.4.2"
        String pipelineUtilModule = "pipeline-util/0.8.20-5-ga0a29bb-foss-2024a"
        String lofreqModule = "LoFreq/2.1.5-foss-2024a"
        Reference reference
        IndexedFile dbsnp
        #IndexedFile cosmic
        IndexedFile gnomadOnlyAfVcf
        IndexedFile panelOfNormalsVcf
        Array[IndexedFile] knownSites
        File? funcotatorDsTar
        #Array[SampleDescriptor] sample
        File? sampleJson
        SampleConfig? sampleConfigIn 
        Array[File] targetIntervals
    }
    #This is a toggle to either use the sampleconfig in the input json or as separate file for conveniance 
    SampleConfig sampleConfig = if (defined(sampleConfigIn)) then select_first([sampleConfigIn]) else read_json(select_first([sampleJson]))
    
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
    
    call somaticLofreq.BamsToSomaticLofreqVariants as somaticLofreq {
        input:
            picardModule = picardModule,
            gatkModule = gatkModule,
            samtoolsModule = samtoolsModule,
            bcftoolsModule = bcftoolsModule,
            pipelineUtilModule = pipelineUtilModule,
            lofreqModule = lofreqModule,
            reference=reference,
            sampleConfigIn=sampleConfig,
            targetIntervals=targetIntervals
    }
    call somaticFreebayes.BamsToSomaticFreebayesVariants as somaticFreebayes {
        input:
            bcftoolsModule = bcftoolsModule,
            freebayesModule = freebayesModule, 
            pipelineUtilModule = pipelineUtilModule,
            reference=reference,
            sampleConfigIn=sampleConfig,
            targetIntervals=targetIntervals
    }
    call somaticMutect.BamsToSomaticMutectVariants as somaticMutect {
        input:
            picardModule = picardModule,
            gatkModule = gatkModule,
            bcftoolsModule = bcftoolsModule,
            pipelineUtilModule = pipelineUtilModule,
            reference = reference,
            dbsnp = dbsnp,
            gnomadOnlyAfVcf = gnomadOnlyAfVcf,
            panelOfNormalsVcf = panelOfNormalsVcf,
            knownSites = knownSites,
            sampleConfigIn = sampleConfig, 
            targetIntervals = targetIntervals
    }
    

    scatter (scatteredtargetsIdx in range(length(targetIntervals))) {
        #Mutect downstream
        call bcftools.Isec as mergeVariantLists {
            input:
                bcftoolsModule=bcftoolsModule,
                inputVcfs=select_all([somaticMutect.vcf[scatteredtargetsIdx],somaticFreebayes.vcf[scatteredtargetsIdx],somaticLofreq.vcf[scatteredtargetsIdx]]),
                inputIdxVcfs=select_all([somaticMutect.idxVcf[scatteredtargetsIdx],somaticFreebayes.idxVcf[scatteredtargetsIdx],somaticLofreq.idxVcf[scatteredtargetsIdx]]),
                outputBasename='variants_merged_scat'+ scatteredtargetsIdx,
        }
        call gatk.HaplotypeCallerRecall as annotateMergedLists {
            input:
                gatkModule=gatkModule,
                reference=reference,
                inputBams=gatherBams.link,
                inputBais=gatherBais.link,
                recallIndexedVcf=mergeVariantLists.idxVcf,
                outputVcfBasename="variants_haplotypecaller_anno_merged_scat"+scatteredtargetsIdx,
                targetScatter = length(targetIntervals)
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
                inputVcfsFiles=[somaticMutect.vcf[scatteredtargetsIdx],somaticFreebayes.vcf[scatteredtargetsIdx],somaticLofreq.vcf[scatteredtargetsIdx]],
                inputVcfs=[somaticMutect.idxVcf[scatteredtargetsIdx],somaticFreebayes.idxVcf[scatteredtargetsIdx],somaticLofreq.idxVcf[scatteredtargetsIdx]],
                outputBasename='annotCallers_merged_scat'+ scatteredtargetsIdx,
        }
        if(defined(funcotatorDsTar)){
            call gatk.Funcotator as funcotateSomaticScattered {
                input:
                    inputVcf=annotateCallers.vcfOut,
                    inputVcfFile=annotateCallers.vcf,
                    reference=reference,
                    funcotatorDsTar=select_first([funcotatorDsTar,reference.fasta]),
                    outputVcfBasename="somatic_funcotated_scat"+ scatteredtargetsIdx,
            }
        }
        File somaticScatteredVcf = select_first([funcotateSomaticScattered.vcf,annotateCallers.vcf])
        IndexedFile somaticScatteredIdxVcf = select_first([funcotateSomaticScattered.vcfOut,annotateCallers.vcfOut])

        #this function below is a band aid due to the level above the scatter not seeing the output as an array.
        call common.CreateIndexedLink as linkSomaticScattered {
            input:
                inputFile=somaticScatteredVcf,
                indexFile=somaticScatteredIdxVcf.index,
                hardLink=true
        }
    }
    
    #merge all regions
    call bcftools.Concat as mergeScatteredRegions {
        input:
            inputVcfs=linkSomaticScattered.file,
            inputIndexedVcfs=linkSomaticScattered.out,
            outputBasename="project_somatic",
            bcftoolsModule=bcftoolsModule,
    }

    
    call picard.SortVcfsIndexed as gatherMutect {
        input:
            picardModule = picardModule,
            inputVcfs=somaticMutect.vcf,
            outputPrefix="project_mutect2",
            createIndex=true
    }
    call picard.SortVcfsIndexed as gatherFreebayes {
        input:
            picardModule = picardModule,
            inputVcfs=somaticFreebayes.vcf,
            outputPrefix="project_freebayes",
            createIndex=true
    }
    call picard.SortVcfsIndexed as gatherLofreq {
        input:
            picardModule = picardModule,
            inputVcfs=somaticLofreq.vcf,
            outputPrefix="project_lofreq",
            createIndex=true
    }

    output {
        #output files of workflow
        #needs to be debugged probably
        Array[File] somaticScatVcf = linkSomaticScattered.file
        Array[IndexedFile] somaticScatIdxVcf = linkSomaticScattered.out
        IndexedFile somaticIdxVcf = mergeScatteredRegions.vcfOut
        #Array[File] bcftoolsStats = [freebayesRawStats.stats,freebayesFilteredStats.stats,mutectFilteredStats.stats,mutectRawStats.stats,lofreqStats.stats]
    }

    meta {
        author: "MMTerpstra"
        description: "This is the single samples bams to variants workflow."
    }  
}