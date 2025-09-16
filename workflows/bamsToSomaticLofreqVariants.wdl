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

workflow BamsToSomaticLofreqVariants {
    input {
        String picardModule = "picard/2.26.10-Java-8-LTS"
        String gatkModule = "GATK/4.2.4.1-Java-8-LTS"
        String samtoolsModule = "SAMtools/1.15.1-GCC-11.3.0"
        String bcftoolsModule = "BCFtools/1.21-GCC-12.2.0"
        #String freebayesModule = "freebayes/1.3.7-gfbf-2024a-R-4.4.2"
        String pipelineUtilModule = "pipeline-util/0.8.20-5-ga0a29bb-foss-2024a"
        String lofreqModule = "LoFreq/2.1.5-foss-2024a"
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
    #this should be the best way of selecting normal samples 
    #idk how it goes for wdl validation though
    #Array[SampleDescriptor] normalSamples = []
    scatter (sample in sampleConfig.samples) {
        call lofreq.LoFreqViterbiCached as viterbiRealingment {
            input:
                inputBamIndexed=select_first([sample.alignedReads]),
                reference=reference,
                picardModule=picardModule,
                lofreqModule=lofreqModule,
                outputBasename=sample.name + '_viterbi'
        }
        
        call common.AddAlignedReadsToSampleDescriptor as addViterbiBamToSample {
            input:
                sample=sample,
                bam=viterbiRealingment.bamOut
        }
        
        SampleDescriptor viterbiSample = addViterbiBamToSample.sampleUpdated
    }

    #this works under the assumtion that at least one sampledescriptor is present for each normal name
    
    scatter(normalname in select_all(uniqueNormals.outArray)){
        scatter (sample in viterbiSample) {
            if (sample.name == normalname){
                SampleDescriptor viterbiNormalSampleScat = sample
            }
        }
        Array[SampleDescriptor] viterbiNormal = select_all(viterbiNormalSampleScat)
    }
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


    Array [SampleDescriptor] viterbiNormalSamples = select_all(flatten(viterbiNormal))
    #this might also be a way of solving the upper problem of merging the optional normalMatched samples to a list of normal samples
    #call common.SelectNormals as selectNormalSamples {
    #    input:
    #        normaleNames = select_all(uniqueNormals.outArray),
    #        samples = sampleConfig.samples
    #}

    scatter (scatteredtargetsIdx in range(length(targetIntervals))) {
        scatter(sample in viterbiSample){
            IndexedFile sampleLofreqIndexedBam = select_first([sample.alignedReads])
            scatter(normal in viterbiNormalSamples){
                IndexedFile normalIndexedBam = select_first([normal.alignedReads])
                if (normal.name != sample.name ){
                    call lofreq.LoFreqSomatic as loFreqSomaticScat {
                        input:
                            lofreqModule=lofreqModule,
                            pipelineUtilModule=pipelineUtilModule,
                            picardModule=picardModule,
                            inputNormalBam=normalIndexedBam.file,
                            inputNormalBamIndex=normalIndexedBam.index,
                            tumorSampleName=sample.name,
                            normalSampleName=normal.name,
                            inputTumorBam=sampleLofreqIndexedBam.file,
                            inputTumorBamIndex=sampleLofreqIndexedBam.index,
                            reference=reference,
                            targetIntervalList=targetIntervals[scatteredtargetsIdx],
                            outputVcfBasename='lofreq_sample_'+sample.name+'_norm_'+normal.name+'_scat'+scatteredtargetsIdx+'_',
                            targetScatter = length(targetIntervals)
                    }
                    call util.Callerise as lofreqSampleCalleriseScat {
                        input:
                            pipelineUtilModule=pipelineUtilModule,
                            inputVcf = loFreqSomaticScat.vcf,
                            outputBase = 'lofreq_caller_sample_'+sample.name+'_norm_'+normal.name+'_scat'+scatteredtargetsIdx,
                            caller="LoFreq_"
                    }
                    call bcftools.Norm as normSampleLofreqScat {
                        input:
                            bcftoolsModule=bcftoolsModule,
                            reference=reference,
                            inputVcf=lofreqSampleCalleriseScat.vcf,
                            outputBasename ='lofreq_norm_sample_'+sample.name+'_norm_'+normal.name+'_scat'+scatteredtargetsIdx,
                    }
                    
                    
                }
            }
            
            Array[File] normalisedTumorNormalSampleLofreqVcfs = select_all(normSampleLofreqScat.vcf)
            Array[IndexedFile] normalisedTumorNormalSampleLofreqIdxVcfs = select_all(normSampleLofreqScat.vcfOut)
            #merge all tumor normals into one for each tumor
            call bcftools.Isec as mergeLoFreqNormalsScat {
                input:
                    bcftoolsModule=bcftoolsModule,
                    inputVcfs=normalisedTumorNormalSampleLofreqVcfs,
                    inputIdxVcfs=normalisedTumorNormalSampleLofreqIdxVcfs,
                    outputBasename='Lofreq_variants_merged_normalisec_scat'+ scatteredtargetsIdx,
                    minIntersect=length(select_all(normSampleLofreqScat.vcf)),
            }
            

            #call gatk.HaplotypeCallerRecall as annotateLofreqNormalMergedLists {
            #    input:
            #        gatkModule=gatkModule,
            #        reference=reference,
            #        inputBams=gatherBams.link,
            #        inputBais=gatherBais.link,
            #        recallIndexedVcf=mergeLoFreqNormalsScat.vcfOut,
            #        outputVcfBasename="lofreq_variants_haplotypecaller_anno_norm_merged_scat"+scatteredtargetsIdx,
            #        targetScatter =  length(targetIntervals)
            #}
            #call bcftools.Norm as normaliseLofreqNormalsMerged {
            #    input:
            #        bcftoolsModule=bcftoolsModule,
            #        reference=reference,
            #        inputVcf=annotateLofreqNormalMergedLists.vcfOut.file,
            #        outputBasename ="norm_merged_scat_"+scatteredtargetsIdx,
            #}
            #this could be better...
            #call util.ReannotateVariants as annotateLofreqNormalsMerged {
            #    input:
            #        pipelineUtilModule=pipelineUtilModule,
            #        combinedVariants=normaliseLofreqNormalsMerged.vcfOut,
            #        inputVcfsFiles=select_all(normSampleLofreqScat.vcf),
            #        inputVcfs=select_all(normSampleLofreqScat.vcfOut),
            #        outputBasename='lofreq_annot_merged_normals_scat'+ scatteredtargetsIdx,
            #}
        }
        Array[File] lofreqScatteredVcfs = select_all(flatten(normalisedTumorNormalSampleLofreqVcfs))
        Array[IndexedFile] lofreqScatteredIdxVcfs = select_all(flatten(normalisedTumorNormalSampleLofreqIdxVcfs))
        #merge all samples into one
        call bcftools.Isec as mergeLoFreqVariants  {
            input:
                bcftoolsModule=bcftoolsModule,
                inputVcfs=select_all(mergeLoFreqNormalsScat.vcf),
                inputIdxVcfs=select_all(mergeLoFreqNormalsScat.vcfOut),
                outputBasename='Lofreq_variants_merged_scat'+ scatteredtargetsIdx,
                minIntersect=1,
        }
        
        call gatk.HaplotypeCallerRecall as annotateLofreqMergedLists {
            input:
                gatkModule=gatkModule,
                reference=reference,
                inputBams=gatherBams.link,
                inputBais=gatherBais.link,
                recallIndexedVcf=mergeLoFreqVariants.vcfOut,
                outputVcfBasename="lofreq_variants_haplotypecaller_anno_merged_scat"+scatteredtargetsIdx,
                targetScatter =  length(targetIntervals)

        }
        call bcftools.Norm as normaliseLofreqMerged {
            input:
                bcftoolsModule=bcftoolsModule,
                reference=reference,
                inputVcf=annotateLofreqMergedLists.vcf,
                outputBasename ="norm_merged_scat_"+scatteredtargetsIdx,
        }
        #this could be better...
        call util.ReannotateVariants as annotateLofreqCallers {
            input:
                pipelineUtilModule=pipelineUtilModule,
                combinedVariants=normaliseLofreqMerged.vcfOut,
                inputVcfsFiles=flatten([lofreqScatteredVcfs,[annotateLofreqMergedLists.vcf]]),
                inputVcfs=flatten([lofreqScatteredIdxVcfs,[annotateLofreqMergedLists.vcfOut]]),
                outputBasename='lofreq_annot_merged_scat'+ scatteredtargetsIdx,
        }
        #filtering for normal samples
        call bcftools.ViewSamples as lofreqNormals {
            input:
                bcftoolsModule=bcftoolsModule,
                inputVcf = annotateLofreqCallers.vcf,
                outputBasename = "freebayes_normals_scat_"+scatteredtargetsIdx,
                samples=uniqueNormals.outArray
        }
        call util.AdFilter as lofreqFilterGermline {
            input:
                inputVariantsToFilter = annotateLofreqCallers.vcfOut,
                inputVcfsFiles=[lofreqNormals.vcf],
                inputVcfs=[lofreqNormals.vcfOut],
                outputBasename='freebayes_normals_filt_scat_'+ scatteredtargetsIdx,
        }
    }
    #calcing lofreq raw values pre merge would be nasty. 
    call bcftools.Stats as lofreqStats {
        input:
            inputVcfs=lofreqFilterGermline.vcf,
            inputIndexedVcfs=lofreqFilterGermline.vcfOut,
            outputBasename="lofreq_stats",
            bcftoolsModule=bcftoolsModule,
    }
    
    output {
        #output files of workflow
        #needs to be debugged probably
        Array[File] vcf = lofreqFilterGermline.vcf
        Array[IndexedFile] idxVcf = lofreqFilterGermline.vcfOut
        Array[File] bcftoolsStats = [lofreqStats.stats]
    }

    meta {
        author: "MMTerpstra"
        description: "This is the single samples bams to variants workflow."
    }  
}