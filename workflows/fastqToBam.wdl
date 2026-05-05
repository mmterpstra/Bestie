version 1.0

import "../structs.wdl" as structs
import "../tasks/common.wdl" as common
import "../tasks/trimgalore.wdl" as trimgalore
import "../tasks/cutadapt.wdl" as cutadapt
import "../tasks/fastqc.wdl" as fastqc
import "../tasks/picard.wdl" as picard
import "../tasks/samtools.wdl" as samtools
import "../tasks/fgbio.wdl" as fgbio
import "../tasks/fgbio_picard.wdl" as fgbio_picard
import "../tasks/alignment.wdl" as align
import "../tasks/gatk.wdl" as gatk
import "../tasks/ichorcna.wdl" as ichorcna
import "../tasks/marktrimming.wdl" as marktrim
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
        String marktrimmingModule = "marktrimming/0.0.2-GCC-12.2.0"
        Boolean runCutadapt = false 
        Boolean removeDuplicates = false
        String cutadaptModule = "cutadapt/4.2-GCCcore-11.3.0"
        Array[String] read1Adapters = ["AGATCGGAAGAGC"]
        Array[String] read2Adapters = ["AGATCGGAAGAGC"]
        Boolean runTwistUmi = false
        Boolean mergeBamFilesCoordinateSort = true
        Boolean runDuplexConsensus = false
        Boolean runBaseQualityRecalibration = true
        Reference reference
        BwaIndex referenceBwaIndex
        IndexedFile dbsnp
        Array[IndexedFile] knownSites
        SampleDescriptor sample
        File targetIntervalList
    }
    
    #Link the sample specific value back to this sample or use the dafault value (usually false)
    Boolean runTwistUmiSample = if (defined(sample.runTwistUmi)) then select_first([sample.runTwistUmi,runTwistUmi]) else runTwistUmi
    Boolean runCutadaptSample = if (defined(sample.runTwistUmi)) then select_first([sample.runTwistUmi,runTwistUmi]) else runCutadapt

    Boolean coordinateSort = if (runTwistUmiSample) then true else mergeBamFilesCoordinateSort
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
        if (defined(rg.fastqUmi)) {
            call common.CreateLink as getfastqUmi {
                input:
                    inputFile = select_first([rg.fastqUmi]),
                    outputPath = sample.name + "_" + rg.flowcell + "_" + rg.identifier + "_umi.fastq.gz"
            }
        }
        call fastqc.FastQCPaired as fastqc {
            input:
                fastqcModule = fastqcModule,
                inputFastq = getfastq1.link,
                inputFastq2 = getfastq2.link,
                outputFastqcBasename = sample.name + "_" + rg.flowcell + "_" + rg.identifier
        }
        
        #there should be more logic to catchall here but atm idk
        Boolean extractUmisFromReadNames = if (defined(rg.extractUmisFromReadNames)) then select_first([rg.extractUmisFromReadNames,false]) else false
        
        call fgbio_picard.FastqToFgUnmappedBamPicardSortedScattered as fastqToUbams {
            input:
                inputFastq1 = getfastq1.link,
                inputFastq2 = getfastq2.link,
                fgbioModule = fgbioModule,
                picardModule = picardModule,
                sampleName = sample.name, 
                platform = rg.platform,
                platformUnit = rg.run + "_" + rg.barcode1 + "+" + select_first([rg.barcode2,'AAAAAA']) + "." + rg.lane,
                library = select_first([rg.library,sample.name + rg.barcode1 + "+" + select_first([rg.barcode2,'AAAAAA']) ]),
                readgroup = rg.run + "_" + rg.flowcell  + "_" + rg.barcode1 + "+" + select_first([rg.barcode2,'AAAAAA']) + "." + rg.lane,
                outputBamBaseDir = rg.run  + "_" + rg.flowcell  + "_" + rg.barcode1 + "+" + select_first([rg.barcode2,'AAAAAA']) + "." + rg.lane,
                inputUmiFastq1 = getfastqUmi.link,
                readStructureFastq1 = rg.readStructureFastq1,
                readStructureFastq2 = rg.readStructureFastq2,
                readStructureFastqUmi = rg.readStructureFastqUmi,
                extractUmisFromReadNames = rg.extractUmisFromReadNames,
        }
        #
        scatter (scatteredUbamsIdx in range(length(fastqToUbams.ubams))) {
            #dump sorted bam reads for trimming for alternate cutadapt workflow
            #outputs ubamToSortedFastq.fastq1gz and select_first(ubamToSortedFastq.fastq2gz)
            File ubam = fastqToUbams.ubams[scatteredUbamsIdx]
            
            if(! runCutadaptSample) {
                call picard.SamToFastq as ubamToSortedFastq {
                    input:
                        inputBam = ubam,
                        picardModule = picardModule,
                        outputFastqDirBase = 'shard' + scatteredUbamsIdx + "run" +rg.run + "_" + rg.flowcell  + "_" + rg.barcode1 + "+" + select_first([rg.barcode2,'AAAAAA']) + "." + rg.lane 
                }
                call trimgalore.TrimGalore as adaptertrim {
                    input:
                        inputFastq1 = ubamToSortedFastq.fastq1gz,
                        inputFastq2 = select_first([ubamToSortedFastq.fastq2gz]),
                        minimumLength = 0,
                        outputFastq1 = 'shard' + scatteredUbamsIdx + 'sample' + sample.name + "_" + rg.flowcell + "_L" + rg.lane + "_" + rg.identifier + "_trim_R1.fastq.gz",
                        outputFastq2 = 'shard' + scatteredUbamsIdx + 'sample' + sample.name + "_" + rg.flowcell + "_L" + rg.lane + "_" + rg.identifier + "_trim_R2.fastq.gz",
                        memoryGb = 1,
                        trimgaloreModule = trimgaloreModule
                }
            }
            if(runCutadaptSample) {
                call cutadapt.CutadaptUbam as cutadaptPe {
                    input:
                        cutadaptModule = cutadaptModule,
                        picardModule = picardModule,
                        samtoolsModule = samtoolsModule,
                        minimumLength = 0,
                        inputUbam = ubam,
                        outputFastq1 = 'shard' + scatteredUbamsIdx + 'sample' + sample.name + "_" + rg.flowcell + "_L" + rg.lane + "_" + rg.identifier + "_cutadapt_R1.fastq.gz",
                        outputFastq2 = 'shard' + scatteredUbamsIdx + 'sample' + sample.name + "_" + rg.flowcell + "_L" + rg.lane + "_" + rg.identifier + "_cutadapt_R2.fastq.gz",
                        read1Adapters = read1Adapters,
                        read2Adapters = read2Adapters
                }                
            }

            File trimmedFastq1=if(runCutadaptSample) then select_first([cutadaptPe.fastq1,adaptertrim.fastq1]) else select_first([adaptertrim.fastq1,cutadaptPe.fastq1])
            File trimmedFastq2=if(runCutadaptSample) then select_first([cutadaptPe.fastq2,adaptertrim.fastq2]) else select_first([adaptertrim.fastq2,cutadaptPe.fastq2])

            #call marktrim.MarkTrimming as markTrimming {
            #    input:
            #        marktrimmingModule = marktrimmingModule,
            #        inputUbam = ubam,
            #        cutadaptFastq1 = select_first([cutadaptPe.fastq1,adaptertrim.fastq1]),
            #        cutadaptFastq2 = select_first([cutadaptPe.fastq2,adaptertrim.fastq2]),
            #        outputTrimmedBamBase = sample.name + "_cutadapt_" + rg.flowcell + "_" + rg.identifier
            #}    
            #if (runTwistUmiSample) {
            #    call fgbio.ExtractUmisFromBam as ExtractUmis {
            #        input:
            #            fgbioModule=fgbioModule,
            #            inputBam=fgbioFastqToUnmappedBam.unalignedBam,
            #            outputBamBasename=rg.run + "_" + rg.flowcell  + "_" + rg.barcode1 + "+" + select_first([rg.barcode2,'AAAAAA']) + "." + rg.lane + "_unaligned_umi"
            #    }
            #}
            #note: consider adding a step to make the trimgalore compatible with the best practices MarkIlluminaAdapters workflow. Though some (old) software just expect the adapters to be removed and not marked.  

            ##map with bwa
            
            call align.bwaMarktrimmingAlignBamSamtoolsCompression as bwaAlignment {
                input:
                    inputUnalignedBam = ubam,
                    cutadaptFastq1 = select_first([cutadaptPe.fastq1,adaptertrim.fastq1]),
                    cutadaptFastq2 = select_first([cutadaptPe.fastq2,adaptertrim.fastq2]),
                    referenceBwaIndex = referenceBwaIndex,
                    reference = reference,
                    bwaModule = bwaModule,
                    picardModule = picardModule,
                    samtoolsModule = samtoolsModule,
                    marktrimmingModule = marktrimmingModule,
                    outputBamBasename = 'shard' + scatteredUbamsIdx + 'sample' + sample.name + "_" + rg.run + "_" + rg.flowcell  + "_" + rg.barcode1 + "+" + select_first([rg.barcode2,'NNNNNNN']) + "." + rg.lane + "_aligned",
                    coordinateSort = coordinateSort
            }
            #call samtools.RemoveNumis as removeNumis {
            #    input:
            #        inputBam = bwaAlignment.bam,
            #        samtoolsModule = samtoolsModule,
            #        outputBasename = 'shard' + scatteredUbamsIdx + 'sample' + sample.name + "_" + rg.run + "_" + rg.flowcell  + "_" + rg.barcode1 + "+" + select_first([rg.barcode2,'NNNNNNN']) + "." + rg.lane + "_aligned_umicleaned"
            #}
        }
        if(runCutadaptSample) {
              
            call fastqc.FastQCSample as fastqcCutadapt_1 {
                input:
                    fastqcModule = fastqcModule,
                    inputFastqGzs = select_all(cutadaptPe.fastq1),
                    outputPrefix = sample.name + "_cutadapt_" + rg.flowcell + "_" + rg.identifier +"_R1"
            }
            if (defined(rg.fastq2)) {
                call fastqc.FastQCSample as fastqcCutadapt_2 {
                    input:
                        fastqcModule = fastqcModule,
                        inputFastqGzs = select_all(cutadaptPe.fastq2),
                        outputPrefix = sample.name + "_cutadapt_" + rg.flowcell + "_" + rg.identifier +"_R2"
                }
            }
        }
    }
    call fastqc.FastQCSample as fastqcSampleR1 {
        input:
            fastqcModule = fastqcModule,
            inputFastqGzs = select_all(getfastq1.link),
            outputPrefix = sample.name + "_R1",
    }
    if (defined(getfastq2.link)) {
        call fastqc.FastQCSample as fastqcSampleR2 {
            input:
                fastqcModule = fastqcModule,
                inputFastqGzs = select_all(getfastq2.link),
                outputPrefix = sample.name + "_R2",
        }
    }
    #remove pcr duplicates // optical
    
    if(runTwistUmiSample){
        String barcodeTag = "RX"
        Boolean duplexUMI = true
    }
    call picard.SortedMarkDuplicates as sortedMarkDups {
        input:
            barcodeTag = barcodeTag,
            duplexUMI = duplexUMI,
            picardModule = picardModule,
            removeDuplicates = removeDuplicates,
            inputBams = flatten(bwaAlignment.bam),
            outputBamBasename = sample.name + '_sort_markdup_umi',
            outputMetrics = sample.name + '.sort_markdup_umi_metrics'
    }
    #
    call qc.bamQualityControl as bamUmiQualityControl {
        input:
            gatkModule = gatkModule,
            picardModule = picardModule,
            fgbioModule = fgbioModule,
            samtoolsModule = samtoolsModule,
            reference = reference,
            inputBam = sortedMarkDups.bam,
            inputBai = select_first([sortedMarkDups.bai]),
            outputPrefix =  sample.name + '_markdup_umi_sort_qc',
            targetIntervalList = targetIntervalList,
            commonVariants = select_first(knownSites),
            byReadGroup = true
    }
    

    File duplicateMarkedBam = sortedMarkDups.bam
    File duplicateMarkedBai = sortedMarkDups.bai
    #runs Duplexconsensus Pipeline

    if(runTwistUmiSample && runDuplexConsensus){
        call picard.MergeSamFiles as mergeBySample{
            input:
                picardModule = picardModule,
                inputBams = flatten(bwaAlignment.bam),
                outputBamBasename = sample.name + '_samplemerged',
        }
        call picard.SortSam as sortMergedSampleBam {
            input: 
                picardModule = picardModule,
                inputBam = mergeBySample.bam,
                outputBamBasename = sample.name + '_notduplicatemarked_sorted',
        }
        call qc.bamQualityControl as bamQualityControlUnMarked {
            input:
            gatkModule = gatkModule,
            picardModule = picardModule,
            fgbioModule = fgbioModule,
            samtoolsModule = samtoolsModule,
            reference = reference,
            inputBam = sortMergedSampleBam.bam,
            inputBai = select_first([sortMergedSampleBam.bai]),
            outputPrefix =  sample.name + '_notduplicatemarked_qc',
            targetIntervalList = targetIntervalList,
            commonVariants = select_first(knownSites),
            byReadGroup = false
        }
        
        call fgbio.GroupReadsByUmi as groupReadsByUmi {
            input:
                fgbioModule = fgbioModule,
                inputBam = mergeBySample.bam,
                outputBamBasename = sample.name + '_before_umi',
        }
        #call fgbio.CallConsensusReads as callConsensusReads {
        #    input:
        #        fgbioModule = fgbioModule,
        #        inputBam = groupReadsByUmi.bam,
        #        outputBamBasename = sample.name + '_consensus_called',
        #}
        call fgbio.CallDuplexConsensusReads as callDuplexConsensusReads {
            input:
                fgbioModule = fgbioModule,
                inputBam = groupReadsByUmi.bam,
                outputBamBasename = sample.name + '_duplex_called',
        }
        
        call picard.SortSam as sortDuplexBam {
        input: 
            picardModule = picardModule,
            inputBam = callDuplexConsensusReads.bam,
            outputBamBasename = sample.name + '_duplex_sorted',
            sortOrder = "queryname"
        }
        
        #optional filterconsensusreads

        call picard.SamToFastq as ubamToSortedDuplexFastq {
            input:
                inputBam = sortDuplexBam.bam,
                picardModule = picardModule,
                outputFastqDirBase = sample.name + '_duplex_sorted_',
        }
        
        
        if(runCutadaptSample) {
            call cutadapt.Cutadapt as cutadaptDuplexPe {
                input:
                    cutadaptModule = cutadaptModule,
                    minimumLength = 0,
                    inputFastq1 = ubamToSortedDuplexFastq.fastq1gz,
                    outputFastq1 = sample.name + "_duplex_cutadapt_R1.fastq.gz",
                    inputFastq2 = select_first([ubamToSortedDuplexFastq.fastq2gz]),
                    outputFastq2 = sample.name + "_duplex_cutadapt_R2.fastq.gz",
                    read1Adapters = read1Adapters,
                    read2Adapters = read2Adapters
            }                
        }
        call align.bwaMarktrimmingAlignBamSamtoolsCompression as bwaDuplexConsensusAlignment {
            input:
                inputUnalignedBam = sortDuplexBam.bam,
                referenceBwaIndex = referenceBwaIndex,
                cutadaptFastq1 = select_first([cutadaptDuplexPe.fastq1,ubamToSortedDuplexFastq.fastq1gz]),
                cutadaptFastq2 = select_first([cutadaptDuplexPe.fastq2,ubamToSortedDuplexFastq.fastq2gz]),
                reference = reference,
                bwaModule = bwaModule,
                picardModule = picardModule,
                marktrimmingModule = marktrimmingModule,
                outputBamBasename = sample.name + "_duplex_aligned",
                coordinateSort = coordinateSort,
                timeMinutes = 20 + ceil(size(sortDuplexBam.bam, "G")) * 120 * 3, #Due to sorting/extra tags (increase in filesize) in the speed decreases a lot.
                umiTags = runTwistUmi
        }
    }
    
    
    #optional indelrealignment
    #wip or skip
    
    #optional basequality score recalibration

    File prebqsrBam = if(runTwistUmiSample && runDuplexConsensus) then select_first([bwaDuplexConsensusAlignment.bam,duplicateMarkedBam]) else duplicateMarkedBam
    File prebqsrBai = if(runTwistUmiSample && runDuplexConsensus) then select_first([bwaDuplexConsensusAlignment.bai,duplicateMarkedBai]) else duplicateMarkedBai

    call qc.bamQualityControl as bamQualityControlPreBqsr {
        input:
        gatkModule = gatkModule,
        picardModule = picardModule,
        fgbioModule = fgbioModule,
        samtoolsModule = samtoolsModule,
        reference = reference,
        inputBam = select_first([prebqsrBam]),
        inputBai = select_first([prebqsrBai]),
        outputPrefix =  sample.name + '_qc',
        targetIntervalList = targetIntervalList,
        commonVariants = select_first(knownSites),
        byReadGroup = false,
        runSamtools = true
    }

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
        call qc.bamQualityControl as recalibratedBamQualityControl {
            input:
                gatkModule = gatkModule,
                picardModule = picardModule,
                fgbioModule = fgbioModule,
                samtoolsModule = samtoolsModule,
                reference = reference,
                inputBam = applyBQSR.bam,
                inputBai = applyBQSR.bai,
                outputPrefix =  sample.name + '_recalibrated_qc',
                targetIntervalList = targetIntervalList,
                commonVariants = select_first(knownSites),
                byReadGroup = true
        }

    }            
    
    File bam = if runBaseQualityRecalibration then select_first([applyBQSR.bam,prebqsrBam]) else prebqsrBam
    File bai = if runBaseQualityRecalibration then select_first([applyBQSR.bai,prebqsrBai]) else prebqsrBai

    output {
        Array [File] fastqcZip = fastqc.zip
        #cutadapt logs
        Array [File?] cutadaptLogs = select_all(flatten([flatten(cutadaptPe.fastq1Log),[cutadaptDuplexPe.fastq1Log]]))
        Array [File?] cutadaptFastqcZip = select_all(flatten([fastqcCutadapt_1.zip,fastqcCutadapt_2.zip]))
        #markduplicates logs
        File markdupLog = sortedMarkDups.metrics
        #bam output
        IndexedFile bam = {
          "file" : bam,
          "index" : bai 
        }
        
        File qcZip = bamQualityControlPreBqsr.qcZip
        File? preUmiQcZip = bamQualityControlUnMarked.qcZip
        File? umiQcZip = bamUmiQualityControl.qcZip
        File? bqsrQcZip = recalibratedBamQualityControl.qcZip
        File? umiFamilySizeHistogram = groupReadsByUmi.familySizeHistogram
        File? basequalityRecalibratonReport = bqsr.recalibrationReport
    }

    meta {
        author: "MMTerpstra"
        description: "This is the single sample fastq to aligned bam workflow."
    }
}