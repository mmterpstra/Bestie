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

workflow FastqToBam {
    input {
        String fastqcModule = "FastQC/0.11.9-Java-11"
        String trimgaloreModule = "Trim_Galore/0.6.6-GCCcore-9.3.0-Python-3.8.2"
        String bwaModule = "BWA/0.7.17-GCCcore-11.3.0"
        String picardModule = "picard/2.26.10-Java-8-LTS"
        String gatkModule = "GATK/4.2.4.1-Java-8-LTS"
        #String hmmcopyutilsModule = "hmmcopy_utils/5911bf69f1-foss-2022a"
        String samtoolsModule = "SAMtools/1.15.1-GCC-11.3.0"
        String fgbioModule = "fgbio/1.3.0"
        Boolean runCutadapt = false 
        String cutadaptModule = "fgbio/1.3.0"
        Array[String] read1Adapters = ["AGATCGGAAGAGC"]
        Array[String] read2Adapters = ["AGATCGGAAGAGC"]
        Boolean runTwistUmi = false
        Boolean runBaseQualityRecalibration = true
        Reference reference
        BwaIndex referenceBwaIndex
        IndexedFile dbsnp
        Array[IndexedFile] knownSites
        SampleDescriptor sample
        File targetIntervalList
    }
    
    #Link the sample specific value back to this sample or use the dafault value (usually false)
    Boolean runTwistUmiSample = if (defined(sample.runTwistUmi)) then select_first([sample.runTwistUmi]) else runTwistUmi

    Boolean coordinateSort = if (runTwistUmiSample) then true else false
    scatter (rg in sample.readgroups) {
        #linking for uniform filenames
        call common.CreateLink as getfastq1 {
            input:
                inputFile = rg.fastq1,
                outputPath = sample.name + "_" + rg.flowcell + "_" + rg.identifier + "_R1.fastq.gz"
        }
        if (defined(rg.fastq2)) {
            call common.CreateLink as getfastq2 {
                input:
                    inputFile = select_first([rg.fastq2]),
                    outputPath = sample.name + "_" + rg.flowcell + "_" + rg.identifier + "_R2.fastq.gz"
            }
        }
        call fastqc.FastQCPaired as fastqc {
            input:
                fastqcModule = fastqcModule,
                inputFastq = getfastq1.link,
                inputFastq2 = getfastq2.link,
                outputFastqcBasename = sample.name + "_" + rg.flowcell + "_" + rg.identifier
        }
        #trim adapters
        call trimgalore.TrimGalore as adaptertrim {
            input:
                inputFastq1 = getfastq1.link,
                inputFastq2 = getfastq2.link,
                outputFastq1 = sample.name + "_" + rg.flowcell + "_L" + rg.lane + "_" + rg.identifier + "_trim_R1.fastq.gz",
                outputFastq2 = sample.name + "_" + rg.flowcell + "_L" + rg.lane + "_" + rg.identifier + "_trim_R2.fastq.gz",
                memoryGb = 1,
                trimgaloreModule = trimgaloreModule
        }
        if(runCutadapt) {
                call cutadapt.Cutadapt as cutadaptPe {
                    input:
                        cutadaptModule = cutadaptModule,
                        inputFastq1 = getfastq1.link,
                        outputFastq1 = sample.name + "_" + rg.flowcell + "_L" + rg.lane + "_" + rg.identifier + "_cutadapt_R1.fastq.gz",
                        inputFastq2 = getfastq2.link,
                        outputFastq2 = sample.name + "_" + rg.flowcell + "_L" + rg.lane + "_" + rg.identifier + "_cutadapt_R2.fastq.gz",
                        read1Adapters = read1Adapters,
                        read2Adapters = read2Adapters
                }
            call fastqc.FastQCPaired as fastqcCutadapt {
                input:
                    fastqcModule = fastqcModule,
                    inputFastq = getfastq1.link,
                    inputFastq2 = getfastq2.link,
                    outputFastqcBasename = sample.name + "_cutadapt_" + rg.flowcell + "_" + rg.identifier
            }
        }        
        #align
        ##to samconversion
        
        call picard.FastqToUnmappedBam as fastqToUnmappedBam {
            input:
                inputFastq1 = select_first([cutadaptPe.fastq1,adaptertrim.fastq1]),
                inputFastq2 = select_first([cutadaptPe.fastq2,adaptertrim.fastq2]),
                picardModule = picardModule,
                sampleName = sample.name, 
                platform = rg.platform,
                platformUnit = rg.run + "_" + rg.barcode1 + "+" + select_first([rg.barcode2,'AAAAAA']) + "." + rg.lane,
                readGroupName = rg.run + "_" + rg.flowcell  + "_" + rg.barcode1 + "+" + select_first([rg.barcode2,'AAAAAA']) + "." + rg.lane,
                outputUnalignedBam = rg.run  + "_" + rg.flowcell  + "_" + rg.barcode1 + "+" + select_first([rg.barcode2,'AAAAAA']) + "." + rg.lane + "_unaligned.bam",
        }
        
        if (runTwistUmiSample) {
            call fgbio.ExtractUmisFromBam as ExtractUmis {
                input:
                    inputBam=fastqToUnmappedBam.unalignedBam,
                    outputBamBasename=rg.run + "_" + rg.flowcell  + "_" + rg.barcode1 + "+" + select_first([rg.barcode2,'AAAAAA']) + "." + rg.lane + "_unaligned_umi.bam",
                    fgbioModule=fgbioModule
            }
        }
        #note: consider adding a step to make the trimgalore compatible with the best practices MarkIlluminaAdapters workflow. Though some (old) software just expect the adapters to be removed and not marked.  

        ##map with bwa

        call align.bwaAlignBam as bwaAlignment {
            input:
                inputUnalignedBam = select_first([ExtractUmis.bam,fastqToUnmappedBam.unalignedBam]),
                referenceBwaIndex = referenceBwaIndex,
                reference = reference,
                bwaModule = bwaModule,
                picardModule = picardModule,
                outputBamBasename = sample.name + rg.run + "_" + rg.flowcell  + "_" + rg.barcode1 + "+" + select_first([rg.barcode2,'NNNNNNN']) + "." + rg.lane + "_aligned",
                coordinateSort = coordinateSort
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

    if(runTwistUmiSample){
        call picard.MergeSamFiles as mergeBySample{
            input:
                picardModule = picardModule,
                inputBams = bwaAlignment.bam,
                outputBamBasename = sample.name + '_sample',
        }
        call fgbio.GroupReadsByUmi as groupReadsByUmi {
            input:
                fgbioModule = fgbioModule,
                inputBam = mergeBySample.bam,
                outputBamBasename = sample.name + '_grouped_by_umi',
        }
        call picard.SortSam as sortMergedSampleBam {
            input: 
                picardModule = picardModule,
                inputBam = mergeBySample.bam,
                outputBamBasename = sample.name + '_sorted',
        }
        call qc.bamQualityControl as bamQualityControlUnMarked {
        #call gatk.CollectMultipleMetrics as CollectMultipleMetrics {
            input:
            gatkModule = gatkModule,
            picardModule = picardModule,
            reference = reference,
            inputBam = sortMergedSampleBam.bam,
            inputBai = sortMergedSampleBam.bai,
            outputPrefix =  sample.name + '_notduplicatemarked_qc',
            targetIntervalList = targetIntervalList,
            byReadGroup = false
        }
        call fgbio.CallDuplexConsensusReads as callDuplexConsensusReads {
            input:
                fgbioModule = fgbioModule,
                inputBam = groupReadsByUmi.bam,
                outputBamBasename = sample.name + '_duplex_called',
        }
        call align.bwaAlignBam as bwaDuplexConsensusAlignment {
            input:
                inputUnalignedBam = callDuplexConsensusReads.bam,
                referenceBwaIndex = referenceBwaIndex,
                reference = reference,
                bwaModule = bwaModule,
                picardModule = picardModule,
                outputBamBasename = sample.name + "_duplex_aligned",
                coordinateSort = coordinateSort
        }

    }
    #remove pcr duplicates
    call picard.MarkDuplicates as markDups {
        input:
            picardModule = picardModule,
            inputBams = bwaAlignment.bam,
            outputBamBasename = sample.name + '_markdup',
            outputMetrics = sample.name + '.markdup_metrics'

    }
    #sort bam by coordinate order
    call picard.SortSam as sortBam {
        input: 
            picardModule = picardModule,
            inputBam = markDups.bam,
            outputBamBasename = sample.name + '_sort'
            
    }
    #optional indelrealignment
    #wip or skip
    
    #optional basequality score recalibration

    File prebqsrBam = if(runTwistUmiSample) then select_first([bwaDuplexConsensusAlignment.bam,sortBam.bam]) else sortBam.bam
    File prebqsrBai = if(runTwistUmiSample) then select_first([bwaDuplexConsensusAlignment.bai,sortBam.bai]) else sortBam.bai
    if(runBaseQualityRecalibration){
        call gatk.BaseQualityScoreRecalibration as bqsr {
            input:
                inputBam=prebqsrBam,
                outputRecalibrationReport=sample.name + '_recalibration.txt',
                reference=reference,
                gatkModule=gatkModule,
                dbsnp=dbsnp,
                knownSites=knownSites
        }
        call gatk.ApplyBQSR as applyBQSR {
            input:
                inputBam=prebqsrBam,
                inputBai=prebqsrBai,
                recalibrationReport=bqsr.recalibrationReport,
                reference=reference,
                gatkModule=gatkModule,
                outputBamBasename=sample.name + '_recalibrated'
        }
    }        
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
            targetIntervalList = targetIntervalList,
            byReadGroup = true
    }
    
    File bam = if runBaseQualityRecalibration then select_first([applyBQSR.bam,prebqsrBam]) else prebqsrBam
    File bai = if runBaseQualityRecalibration then select_first([applyBQSR.bai,prebqsrBai]) else prebqsrBai

    output {
        Array [File] fastqcZip = fastqc.zip
        #cutadapt logs
        Array [File?] cutadaptLogs = select_all(cutadaptPe.fastq1Log)
        #markduplicates logs
        File markdupLog = markDups.metrics
        #bam output
        IndexedFile bam = {
          "file" : bam,
          "index" : bai 
        }
        File qcZip = bamQualityControl.qcZip
        File? umiQcZip = bamQualityControlUnMarked.qcZip
    }

    meta {
        author: "MMTerpstra"
        description: "This is the single samples fastqs to aligned bam workflow."
    }
}
