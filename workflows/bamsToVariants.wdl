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
import "../tasks/vep.wdl" as vep
import "../workflows/qc.wdl" as qc
import "../workflows/bamsToGermlineVariants.wdl" as bamsToGermlineVariants
import "../workflows/bamsToSomaticVariants.wdl" as bamsToSomaticVariants



workflow BamsToVariants {
    input {
        String picardModule = "picard/2.26.10-Java-8-LTS"
        String gatkModule = "GATK/4.2.4.1-Java-8-LTS"
        String samtoolsModule = "SAMtools/1.15.1-GCC-11.3.0"
        String bcftoolsModule = "BCFtools/1.21-GCC-12.2.0"
        String freebayesModule = "freebayes/1.3.7-gfbf-2024a-R-4.4.2"
        String pipelineUtilModule = "pipeline-util/0.8.20-5-ga0a29bb-foss-2024a"
        String lofreqModule = "LoFreq/2.1.5-foss-2024a"
        String vepModule = "VEP/113.3-GCC-13.3.0" 
        Reference reference
        IndexedFile dbsnp
        #IndexedFile cosmic
        IndexedFile gnomadOnlyAfVcf
        IndexedFile panelOfNormalsVcf
        Array[IndexedFile] knownSites
        File? funcotatorDsTar
        File? vepCacheTar
        #Array[SampleDescriptor] sample
        File? sampleJson
        SampleConfig? sampleConfigIn 
        File targetIntervalList
        Int targetScatter
    }
    #This is a toggle to either use the sampleconfig in the input json or as separate file for conveniance 
    SampleConfig sampleConfig = if (defined(sampleConfigIn)) then select_first([sampleConfigIn]) else read_json(select_first([sampleJson]))

    call picard.SplitAndPadIntervals as splitIntervals {
        input:
            picardModule = picardModule,
            inputIntervalListFile = targetIntervalList,
            outputPrefix = "Targets_padded_and_split",
            targetScatter = targetScatter,
            padding = 200
    }
    
    call bamsToSomaticVariants.BamsToSomaticVariants as somaticScattered {
        input:
            picardModule = picardModule,
            gatkModule = gatkModule,
            samtoolsModule = samtoolsModule,
            bcftoolsModule = bcftoolsModule,
            freebayesModule = freebayesModule,
            pipelineUtilModule = pipelineUtilModule,
            lofreqModule = lofreqModule,
            vepModule = vepModule,
            reference = reference,
            dbsnp = dbsnp,
            gnomadOnlyAfVcf = gnomadOnlyAfVcf,
            funcotatorDsTar = funcotatorDsTar,
            vepCacheTar = vepCacheTar,
            panelOfNormalsVcf = panelOfNormalsVcf,
            knownSites = knownSites,
            sampleConfigIn = sampleConfig, 
            targetIntervals = splitIntervals.paddedScatteredIntervalList
    }

    if(defined(vepCacheTar)){
        call vep.AnnotateTmpStorage as vepannotateSomaticScattered {
            input:
                vepModule = vepModule,
                reference = reference,
                vepCacheTar = select_first([vepCacheTar,reference.fasta]),
                inputIdxVcf = somaticScattered.somaticIdxVcf,
                outputBasename = "somatic_vep",
        }

    }
    #call bamsToGermlineVariants.BamsToGermlineVariants as germline {
    #    input:
    #        gatkModule = gatkModule,
    #        reference = reference,
    #        sampleConfigIn = sampleConfig,
    #        targetIntervals = splitIntervals.paddedScatteredIntervalList
    #}
    #Array[IndexedFile] sampleIndexedBams = sampleConfig[].samples.alignedReads

    output {
        #output files of workflow
        #needs to be debugged probably
        IndexedFile somaticIdxVcf = somaticScattered.somaticIdxVcf
        #IndexedFile haplotypecallerGVcf = germline.idxGvcf
        #IndexedFile haplotypecallerIdxVcf = germline.idxVcf
        #Array[File] bcftoolsStats = [freebayesRawStats.stats,freebayesFilteredStats.stats,mutectFilteredStats.stats,mutectRawStats.stats,lofreqStats.stats]
    }

    meta {
        author: "MMTerpstra"
        description: "This is the single samples bams to variants workflow."
    }  
}