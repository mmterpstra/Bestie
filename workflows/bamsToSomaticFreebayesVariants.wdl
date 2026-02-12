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


workflow BamsToSomaticFreebayesVariants {
    input {
        #String picardModule = "picard/2.26.10-Java-8-LTS"
        String gatkModule = "GATK/4.2.4.1-Java-8-LTS"
        #String samtoolsModule = "SAMtools/1.15.1-GCC-11.3.0"
        String bcftoolsModule = "BCFtools/1.21-GCC-12.2.0"
        String freebayesModule = "freebayes/1.3.7-gfbf-2024a-R-4.4.2"
        String pipelineUtilModule = "pipeline-util/0.8.20-5-ga0a29bb-foss-2024a"
        #String lofreqModule = "LoFreq/2.1.5-foss-2024a"
        Reference reference
        #IndexedFile dbsnp
        #IndexedFile cosmic
        #IndexedFile gnomadOnlyAfVcf
        #IndexedFile panelOfNormalsVcf
        #Array[IndexedFile] knownSites
        #File? funcotatorDsTar
        #Array[SampleDescriptor] sample
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
    call common.InverseSelection as uniqueTumors {
        input:
            array = sampleNames,
            selection = uniqueNormals.outArray
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
        
        call freebayes.FreebayesSomatic as freebayesSomatic {
            input:
                freebayesModule = freebayesModule,
                reference = reference,
                inputBams = gatherBams.link,
                inputBamIndexes = gatherBais.link,
                targetIntervalList = targetIntervals[scatteredtargetsIdx],
                outputVcfBasename = "freebayes_joint_calls_scat"+scatteredtargetsIdx,
                targetScatter = length(targetIntervals)
        }
        #fixing sample order
        call bcftools.ViewSamples as fixFreebayesSomatic {
            input:
                bcftoolsModule=bcftoolsModule,
                inputVcf=freebayesSomatic.vcf,
                samples=sampleNames,
                outputBasename ='freebayes_joint_calls_samplenames_fixed_scat'+scatteredtargetsIdx,
        }
    
        #freebayes downstream filter for normal sample callings
        call bcftools.ViewSamples as freebayesNormals {
            input:
                bcftoolsModule=bcftoolsModule,
                inputVcf = fixFreebayesSomatic.vcf,
                outputBasename = "freebayes_normals_scat_"+scatteredtargetsIdx,
                samples=uniqueNormals.outArray
        }

        scatter (tumorSample in uniqueTumors.result){
            call bcftools.ViewSamples as freebayesTumor {
                input:
                    bcftoolsModule=bcftoolsModule,
                    inputVcf = fixFreebayesSomatic.vcf,
                    outputBasename = "freebayes_"+tumorSample+"_scat_"+scatteredtargetsIdx,
                    samples=[tumorSample]
            }
            call util.AdFilter as freebayesTumorFilterGermline {
                input:
                    pipelineUtilModule=pipelineUtilModule,
                    inputVariantsToFilter = freebayesTumor.vcfOut,
                    freq=0.01,
                    normalFreq=0.1,
                    count=5, 
                    inputVcfsFiles=[freebayesNormals.vcf],
                    inputVcfs=[freebayesNormals.vcfOut],
                    outputBasename="freebayes_"+tumorSample+"_normalfilt_scat_"+scatteredtargetsIdx,
            }
        }
        call bcftools.Isec as combineFreebayesTumor {
            input:
                bcftoolsModule=bcftoolsModule,
                inputVcfs=freebayesTumorFilterGermline.vcf,
                inputIdxVcfs=freebayesTumorFilterGermline.vcfOut,
                outputBasename='freebayes_somatic_tumor_iseq_scat'+ scatteredtargetsIdx,
                minIntersect=1,
        }
        #call gatk.FreebayesRecall as recallSomatic {
        #    input:
        #        freebayesModule = freebayesModule,
        #        reference = reference,
        #       inputBams = gatherBams.link,
        #        inputBamIndexes = gatherBais.link,
        #        outputVcfBasename = "freebayes_joint_calls_scat"+scatteredtargetsIdx,
        #        targetScatter = length(targetIntervals),
        #        inputVariants = combineFreebayesTumor.idxVcf,
        #}
        #call bcftools.ViewSamples as fixFreebayesRecallSomatic {
        #    input:
        #        bcftoolsModule=bcftoolsModule,
        #        inputVcf=recallSomatic.vcf,
        #        samples=sampleNames,
        #        outputBasename ='freebayes_recall_samplenames_fixed_scat'+scatteredtargetsIdx,
        #}
        call gatk.HaplotypeCallerRecall as annotateFreebayesMergedLists {
            input:
                gatkModule=gatkModule,
                reference=reference,
                inputBams=gatherBams.link,
                inputBais=gatherBais.link,
                recallIndexedVcf=combineFreebayesTumor.idxVcf,
                outputVcfBasename="freebayes_variants_haplotypecaller_anno_scat"+scatteredtargetsIdx,
                targetScatter =  length(targetIntervals)

        }   
        
        #filterfreebayes?
        call util.Callerise as freebayesScatCallerise {
            input:
                pipelineUtilModule=pipelineUtilModule,
                inputVcf = annotateFreebayesMergedLists.vcf,
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
        #Viewsamples fixes sample order by forcing a sample order to the output samples.
        call bcftools.ViewSamples as fixNormaliseFreeBayes {
            input:
                bcftoolsModule=bcftoolsModule,
                inputVcf=normaliseFreeBayes.vcf,
                samples=sampleNames,
                outputBasename ='freebayes_norm_fixed_scat_'+scatteredtargetsIdx,
        }
    }
    call bcftools.Stats as freebayesRawStats {
        input:
            inputVcfs=fixFreebayesSomatic.vcf,
            inputIndexedVcfs=fixFreebayesSomatic.vcfOut,
            outputBasename="freebayes_raw",
            bcftoolsModule=bcftoolsModule,
    }

    call bcftools.Stats as freebayesFilteredStats {
        input:
            inputVcfs=fixNormaliseFreeBayes.vcf,
            inputIndexedVcfs=fixNormaliseFreeBayes.vcfOut,
            outputBasename="freebayes_filtered",
            bcftoolsModule=bcftoolsModule,
    }

    output {
        #output files of workflow
        #needs to be debugged probably
        Array[File] vcf = fixNormaliseFreeBayes.vcf
        Array[IndexedFile] idxVcf = fixNormaliseFreeBayes.vcfOut
        Array[File] bcftoolsStats = [freebayesRawStats.stats,freebayesFilteredStats.stats]
    }

    meta {
        author: "MMTerpstra"
        description: "This is the single samples bams to somatic variants workflow using freebayes."
    }  
}