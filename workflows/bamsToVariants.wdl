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
import "../tasks/ichorcna.wdl" as ichorcna
import "../workflows/qc.wdl" as qc

workflow BamsToVariants {
    input {
        String picardModule = "picard/2.26.10-Java-8-LTS"
        String gatkModule = "GATK/4.2.4.1-Java-8-LTS"
        String samtoolsModule = "SAMtools/1.15.1-GCC-11.3.0"
        Reference reference
        BwaIndex referenceBwaIndex
        IndexedFile dbsnp
        IndexedFile cosmic
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
    scatter (sample in sampleConfig.samples) {

        scatter (scatteredtargetsIdx in range(length(splitIntervals.paddedScatteredIntervalList))) {
            IndexedFile sampleIndexedBam = select_first([sample.alignedReads])
        
            #haplotypecallergvcf
            call gatk.HaplotypeCallerGVcf as haplotypeCallerGvcf {
                input:
                    gatkModule = gatkModule,
                    reference = reference,
                    inputBam = sampleIndexedBam,
                    outputVcfBasename = sample.name + ".idx_" + scatteredtargetsIdx,
                    targetIntervalList = splitIntervals.paddedScatteredIntervalList[scatteredtargetsIdx]
            }
            
            #mutect
            scatter (sampleControl in sampleConfig.samples) {
                if (defined(sample.control) ){
                    if(sample.control == sampleControl.name){
                        #call gatk.Mutect2 as mutect2 {
                        #    input:
                        #        gatkModule = gatkModule,
                        #        reference = reference,
                        #        inputBam = sampleIndexedBam,
                        #}
                    }
                }
            }
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
    #CombineGVCFs to create a single project gvcf
    call gatk.CombineGVCFs as gatherSamples {
                input:
                    gatkModule = gatkModule,
                    reference = reference,
                    inputGVcfs = gatherRegions.vcfOut,
                    inputGVcfsFiles = gatherRegions.vcf,
                    outputVcfBasename = "haplotypecallergvcf_allsamples",
    }
    #genotype gvcf
    call gatk.GenotypeGVCFs as genotype {
                input:
                    gatkModule = gatkModule,
                    reference = reference,
                    inputGVcfs = gatherRegions.vcfOut,
                    inputGVcfsFiles = gatherRegions.vcf,
                    outputVcfBasename = "haplotypecallergvcf_allsamples",
    }
    
}