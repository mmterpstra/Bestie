version 1.0

import "../structs.wdl" as structs

task CollectMultipleMetrics {
    input {
        File inputBam
        String outputMetricsBasename
        Reference reference
        String gatkModule = "GATK"
        Int? memoryGb = "4"
        Int? javaXmxMemoryMb = floor(memoryGb*0.9*1024)
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
        Int disk = ceil(size(inputBam, "M")*1.1)
        Boolean byReadGroup = false
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
        ml ~{gatkModule}
        gatk --java-options "-Xmx~{javaXmxMemoryMb}m -XX:ParallelGCThreads=1" CollectMultipleMetrics \
        --REFERENCE_SEQUENCE ~{reference.fasta} \
        --INPUT ~{inputBam} \
        --OUTPUT ~{outputMetricsBasename} \
        ~{if byReadGroup then "--METRIC_ACCUMULATION_LEVEL READ_GROUP " else ""} 

    }
    
    output {
        File alignmentMetrics = outputMetricsBasename + ".alignment_summary_metrics"
        File baseDistributionMetrics = outputMetricsBasename + ".base_distribution_by_cycle_metrics"
        File baseDistributionPdf = outputMetricsBasename + ".base_distribution_by_cycle.pdf"
        File insertSizeMetrics = outputMetricsBasename + ".insert_size_metrics"
        File insertSizePdf = outputMetricsBasename + ".insert_size_histogram.pdf"
        File qualityByCycleMetrics = outputMetricsBasename + ".quality_by_cycle_metrics"
        File qualityByCyclePdf = outputMetricsBasename + ".quality_by_cycle.pdf"
        File qualityDistributionMetrics = outputMetricsBasename + ".quality_distribution_metrics"
        File qualityDistributionPdf = outputMetricsBasename + ".quality_distribution.pdf"
        File readLengthPdf = outputMetricsBasename + ".read_length_histogram.pdf"
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}

task CollectHsMetrics {
    input {
        File inputBam
        File targetIntervalList
        String outputMetricsBasename
        Reference reference
        String gatkModule = "GATK"
        Int? memoryGb = "4"
        Int? javaXmxMemoryMb = floor(memoryGb*0.9*1024)
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
        Int disk = ceil(size(inputBam, "M")*2.1)
        Boolean byReadGroup = false
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
        ml ~{gatkModule}

        gatk --java-options "-Xmx~{javaXmxMemoryMb}m -XX:ParallelGCThreads=1" CollectHsMetrics \
        --REFERENCE_SEQUENCE "~{reference.fasta}" \
        --BAIT_INTERVALS "~{targetIntervalList}" \
        --TARGET_INTERVALS "~{targetIntervalList}" \
        --INPUT "~{inputBam}" \
        --OUTPUT "~{outputMetricsBasename}.hs_metrics" \
        --COVERAGE_CAP 10000 \
        ~{if byReadGroup then "--METRIC_ACCUMULATION_LEVEL READ_GROUP " else ""} 

    }
    
    output {
        File hsMetrics = outputMetricsBasename + ".hs_metrics"
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}

task CollectWgsMetrics {
    input {
        File inputBam
        String outputMetricsBasename
        Reference reference
        String gatkModule = "GATK"
        Int? memoryGb = "4"
        Int? javaXmxMemoryMb = floor(memoryGb*0.9*1024)
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
        Int disk = ceil(size(inputBam, "M")*2.1)
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
        ml ~{gatkModule}

        gatk --java-options "-Xmx~{javaXmxMemoryMb}m -XX:ParallelGCThreads=1" CollectWgsMetrics \
        --REFERENCE_SEQUENCE "~{reference.fasta}" \
        --INPUT "~{inputBam}" \
        --OUTPUT "~{outputMetricsBasename}.wgs_metrics" \
        --COVERAGE_CAP 500 
    }
    
    output {
        File wgsMetrics = outputMetricsBasename + ".wgs_metrics"
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}
task DepthOfCoverage {
    input {
        File inputBam
        File targetIntervalList
        String outputMetricsBasename
        Array[Int] summaryCoverageThresholds = [10,20,50,100,200,500,1000,2000,5000,10000]
        Reference reference
        String gatkModule = "GATK"
        Int? memoryGb = "6"
        Int? javaXmxMemoryMb = floor(memoryGb*0.9*1024)
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
        Int disk = ceil(size(inputBam, "M")*2.1)
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
        ml ~{gatkModule}

        gatk --java-options "-Xmx~{javaXmxMemoryMb}m -XX:ParallelGCThreads=1" \
        DepthOfCoverage \
        -R "~{reference.fasta}" \
        --count-type COUNT_READS \
        --output-format TABLE \
        -I ~{inputBam} \
        --output "~{outputMetricsBasename}.dcov_metrics" \
        -L "~{targetIntervalList}" \
        --start 1 \
        --stop 2000 \
        --summary-coverage-threshold ~{sep=' --summary-coverage-threshold ' summaryCoverageThresholds} 
    }
    
    output {
        File dcovMetrics = outputMetricsBasename + ".dcov_metrics"
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}

#wip
task BaseQualityScoreRecalibration {
    input {
        File inputBam
        String outputRecalibrationReport
        String gatkModule = "GATK"
        #File targetIntervalList
        #String outputMetricsBasename
        Reference reference
        IndexedFile dbsnp
        Array [IndexedFile] knownSites
        String gatkModule = "GATK"
        Int? memoryGb = "6"
        Int? javaXmxMemoryMb = floor(memoryGb*0.9*1024)
        Int disk = ceil(size(inputBam, "M")*1.1)
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
        #rray [File] knownSitesVcfs = select_all(knownSites).file
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command <<<
        ml ~{gatkModule}
        KNOWNSITES=$(cat ~{write_objects(knownSites)}| \
            cut -f 1 | \
            tail -n+2)
        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xmx~{javaXmxMemoryMb}m" BaseRecalibrator \
            --reference "~{reference.fasta}" \
            --input "~{inputBam}" \
            --use-original-qualities \
            --output ~{outputRecalibrationReport} \
            --known-sites ~{dbsnp.file} \
            $(printf ' --known-sites %s ' $(printf '%s\n' ${KNOWNSITES[@]}))
            
            #-L Interval list separated format
    >>>
    
    output {
        File recalibrationReport = outputRecalibrationReport
    }

    runtime {
        memory: select_first([memoryGb * 1024,6*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}

#gatherbsqr reports

task ApplyBQSR {
  input {
    File inputBam
    File inputBai
    String gatkModule = "GATK"
    String outputBamBasename
    File recalibrationReport
    Reference reference
    Int? memoryGb = "4"
    Int? javaXmxMemoryMb = floor(memoryGb*0.9*1024)
    Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
    Int disk = ceil(size(inputBam, "M")*2.6)
  }
  #Float referenceSize = size(reference.fasta, "GiB") + size(reference.dict, "GiB") + size(reference.fai, "GiB")
  #Int DiskSize = ceil((size(inputBam, "GiB") * 3 ) + referenceSize) + 20
  command {
    ml ~{gatkModule}
    gatk --java-options "-XX:+PrintFlagsFinal \
        -XX:+PrintGCDetails -Xloggc:gc_log.log \
        -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m -Xmx~{javaXmxMemoryMb}m" \
        ApplyBQSR \
        --create-output-bam-md5 true\
        --create-output-bam-index true\
        --add-output-sam-program-record \
        -R ~{reference.fasta} \
        -I ~{inputBam} \
        --use-original-qualities \
        -O ~{outputBamBasename}.bam \
        -bqsr ~{recalibrationReport} 
  }
  output {
    File bam = outputBamBasename + ".bam"
    File bai = outputBamBasename + ".bai"
    File bammd5sum = outputBamBasename + ".bam.md5"
  }
  runtime {
    memory: select_first([memoryGb * 1024,4*1024])
    timeMinutes: timeMinutes
    disk: disk
  }
}

task HaplotypeCallerGVcf {
    input {
        IndexedFile inputBam
        #File inputBai
        File targetIntervalList
        String outputVcfBasename
        Reference reference
        String gatkModule = "GATK"
        Int memoryGb = "4"
        Int javaMemoryGb = memoryGb - 1
        Boolean makeGvcf = true
        Boolean makeBamOut = true
        Boolean useSpanningEventGenotyping = false
        Float? contamination = 0
        Int targetScatter = 1
        Int timeMinutes = 1 + ceil(size(inputBam.file, "G")) * 120 / targetScatter
        Int? javaXmxMemoryMb = floor(memoryGb*0.9*1024)
        Int disk = ceil(size(inputBam.file, "M")*1.2)
    }

    String vcfSuffix =  if makeGvcf then ".g.vcf.gz" else ".vcf.gz"
    String bamoutArg =  if makeBamOut then "-bamout " + outputVcfBasename + ".bamout.bam" else ""
    File bamIn = inputBam.file
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
        set -e

        ml ~{gatkModule}
        gatk --java-options "-Xmx${javaXmxMemoryMb}m -Xms${javaXmxMemoryMb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            HaplotypeCaller \
            -R ~{reference.fasta} \
            -I ~{bamIn} \
            -L ~{targetIntervalList} \
            -O "~{outputVcfBasename}~{vcfSuffix}" \
            -contamination ~{default=0 contamination} \
            -G StandardAnnotation -G StandardHCAnnotation ~{true="-G AS_StandardAnnotation" false="" makeGvcf} \
            ~{false="--disable-spanning-event-genotyping" true="" useSpanningEventGenotyping} \
            -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
            ~{true="-ERC GVCF" false="" makeGvcf} \
            ~{bamoutArg}

        touch ~{outputVcfBasename}.bamout.bam
    }
    
    output {
        File vcf = outputVcfBasename + vcfSuffix
        File vcfIdx = outputVcfBasename + vcfSuffix + ".tbi"
        IndexedFile vcfOut = {
          "file" : outputVcfBasename + vcfSuffix,
          "index" : outputVcfBasename + vcfSuffix + ".tbi"
        }
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}
# add genomicsdb import#
#--merge-contigs-into-num-partitions 25

task CombineGVCFs {
    input{
        Array[IndexedFile] inputGVcfs
        Array[File] inputGVcfsFiles
        String outputVcfBasename
        String vcfSuffix = ".g.vcf.gz"
        Reference reference
        String gatkModule = "GATK"
        Int memoryGb = "4"
        Int javaXmxMemoryMb = ceil((memoryGb - 0.5) * 1024)
        Int timeMinutes = 1 + ceil(size(inputGVcfsFiles, "G")) * 120
        Int disk = 1 + ceil(size(inputGVcfsFiles, "G") * 2.1) #worst case
      }
    #Array[File] gvcfs = select_all(inputGVcfs)[]["file"]
    command <<<
        ml ~{gatkModule}
         gatk --java-options "-Xmx~{javaXmxMemoryMb}m -Xms~{javaXmxMemoryMb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            CombineGVCFs \
            -R ~{reference.fasta} \
            --variant ~{sep=' --variant ' inputGVcfsFiles} \
            -O ~{outputVcfBasename}~{vcfSuffix}
    >>>
    output {
        File vcf = outputVcfBasename + vcfSuffix
        File vcfIdx = outputVcfBasename + vcfSuffix + ".tbi"
        IndexedFile vcfOut = {
          "file" : outputVcfBasename + vcfSuffix,
          "index" : outputVcfBasename + vcfSuffix + ".tbi"
        }
    }
    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }

}
task GenotypeGVCFs {
    input{
        IndexedFile inputGVcf
        File inputGVcfsFile
        String outputVcfBasename
        String vcfSuffix = ".vcf.gz"
        Reference reference
        String gatkModule = "GATK"
        Int memoryGb = "4"
        Int javaXmxMemoryMb = ceil((memoryGb - 0.5) * 1024)
        Int timeMinutes = 1 + ceil(size(inputGVcfsFile, "G")) * 120
        Int disk = 1 + ceil(size(inputGVcfsFile, "G")*1024 * 1.1) #worst case
    }
    
    #Array[File] gvcfs = select_all(inputGVcfs)[]["file"]
    command <<<
        ml ~{gatkModule}
         gatk --java-options "-Xmx~{javaXmxMemoryMb}m -Xms~{javaXmxMemoryMb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            CombineGVCFs \
            -R ~{reference.fasta} \
            --variant ~{inputGVcfsFile} \
            -O ~{outputVcfBasename}~{vcfSuffix}
    >>>
    output {
        File vcf = outputVcfBasename + vcfSuffix
        File vcfIdx = outputVcfBasename + vcfSuffix + ".tbi"
        IndexedFile vcfOut = {
          "file" : outputVcfBasename + vcfSuffix,
          "index" : outputVcfBasename + vcfSuffix + ".tbi"
        }
    }
    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }

}

task GenomicsDBImport {
    input{
        Array[IndexedFile] inputIdxVcfs
        Array[File] inputVcfs
        String outputBasename
        File interval_list
        Reference reference
        String gatkModule = "GATK"
        Int threads = 6
        Int memoryGb = "5"
        Int javaXmxMemoryMb = ceil((memoryGb - 0.5) * 1024)
        Int timeMinutes = 1 + ceil(size(inputVcfs, "G")) * 120
        Int disk = 1 + ceil(size(inputVcfs, "G") * 2.1) #worst case
      }
    #Array[File] gvcfs = select_all(inputGVcfs)[]["file"]
    command <<<
        set -euo pipefail
        ml ~{gatkModule}
        gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
            --reference ~{reference.fasta} \
            --variant ~{sep=' --variant ' inputVcfs} \
            --genomicsdb-workspace-path ~{outputBasename}.genomicsdb \
            --tmp-dir $TMPDIR \
            --genomicsdb-shared-posixfs-optimizations true \
            --batch-size 50 \
            --reader-threads ~{threads - 1} \
            --merge-contigs-into-num-partitions 25 \
            --merge-input-intervals \
            --intervals ~{interval_list}
        
        tar -cf ~{outputBasename}.genomicsdb.tar ~{outputBasename}.genomicsdb
    >>>
    output {
        File genomicsDbTar = outputBasename + ".genomicsdb.tar" 
    }
    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }

}

task SelectVariants {
    input{
        File genomicsDbTar
        String outputVcfBasename
        String vcfSuffix = ".vcf.gz"
        Reference reference
        String gatkModule = "GATK"
        Int memoryGb = "5"
        Int javaXmxMemoryMb = ceil((memoryGb - 0.5) * 1024)
        Int timeMinutes = 1 + ceil(size(genomicsDbTar, "G")) * 120
        Int disk = 1 + ceil(size(genomicsDbTar, "G") * 2.1) #worst case
    }
    #Array[File] gvcfs = select_all(inputGVcfs)[]["file"]
    command <<<
        set -euo pipefail
        ml ~{gatkModule}
        tar -xf ~{genomicsDbTar}

        gatk --java-options "-Xmx4g -Xms4g" SelectVariants \
            --reference ~{reference.fasta} \
            --variant gendb://"$(basename ~{genomicsDbTar} .tar)" \
            --genomicsdb-workspace-path \
            --tmp-dir=$TMPDIR \
            --output ~{outputVcfBasename}.vcf.gz
    >>>
    output {
        File vcf = outputVcfBasename + ".vcf.gz" 
    }
    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }

}

#
task Funcotator {
    input{
        IndexedFile inputVcf
        File inputVcfFile
        File funcotatorDsTar
        String outputVcfBasename
        String vcfSuffix = ".vcf.gz"
        Reference reference
        String gatkModule = "GATK"
        Int memoryGb = "4"
        Int javaXmxMemoryMb = ceil((memoryGb - 0.5) * 1024)
        Int timeMinutes = 1 + ceil(size(inputVcfFile, "G")) * 120
        Int disk = 1 + ceil(size([inputVcfFile,funcotatorDsTar], "G")*1024 * 1.1) #worst case
    }
    
    #Array[File] gvcfs = select_all(inputGVcfs)[]["file"]
    command <<<
        ml ~{gatkModule}
        mkdir -p $TMPDIR/funcotatorDataSource
        (cd $TMPDIR/funcotatorDataSource && tar -xf ~{funcotatorDsTar})

         gatk --java-options "-Xmx~{javaXmxMemoryMb}m -Xms~{javaXmxMemoryMb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            Funcotator \
                --variant ~{inputVcfFile} \
                --reference  ~{reference.fasta} \
                --ref-version hg38 \
                --data-sources-path $TMPDIR/funcotatorDataSource/* \
                --output ~{outputVcfBasename}~{vcfSuffix} \
                --output-file-format VCF
    >>>
    output {
        File vcf = outputVcfBasename + vcfSuffix
        File vcfIdx = outputVcfBasename + vcfSuffix + ".tbi"
        IndexedFile vcfOut = {
          "file" : outputVcfBasename + vcfSuffix,
          "index" : outputVcfBasename + vcfSuffix + ".tbi"
        }
    }
    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }

}

task HaplotypeCallerRecall {
    input {
        #rray[IndexedFile] inputIndexedBam
        Array[File] inputBams
        Array[File] inputBais
        IndexedFile recallIndexedVcf
        String outputVcfBasename
        Reference reference
        String gatkModule = "GATK"
        Int memoryGb = "12"
        #Int javaMemoryGb = memoryGb - 1
        Boolean makeGvcf = true
        Boolean makeBamOut = false
        #useSpanningEventGenotyping creates spanning event or  aternative reference calls in your dataset usually confusing downstream processing tools 
        Boolean useSpanningEventGenotyping = false
        Float? contamination = 0
        Int targetScatter = 1
        Int timeMinutes = 5 + ceil(size(inputBams, "G")) * 120 / targetScatter
        Int javaXmxMemoryMb = select_first ([floor((memoryGb-1)*0.9*1024),floor(10*0.9*1024)])
        Int disk = ceil(size(inputBams, "M")*1.2)  
    }

    String vcfSuffix =  if makeGvcf then ".g.vcf.gz" else ".vcf.gz"
    String bamoutArg =  if makeBamOut then "-bamout " + outputVcfBasename + ".bamout.bam" else ""
    
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command <<<
        set -e

         >&2 echo " ## "$(date)" ## Localising files" 
        cat ~{write_lines(select_all(inputBams))} \
            ~{write_lines(inputBais)} | \
            (while read FILE; do 
                if [ ! -e "$TMPDIR/""$(basename "$FILE")" ]; then
                    cp "$FILE" "$TMPDIR/"
                    if [[ $FILE =~ \.bam$ ]]; then
                        echo "$TMPDIR/""$(basename "$FILE")">>"./bams_inputs.list" 
                    fi
                else
                    echo "Duplicate file basename spotted $FILE" && exit 1 
                fi
            done )
        >&2 echo " ## "$(date)" ## Loading modules" 

        ml ~{gatkModule}
        gatk --java-options "-Xmx~{javaXmxMemoryMb}m -Xms~{javaXmxMemoryMb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            HaplotypeCaller \
            -R $(realpath ~{reference.fasta}) \
            $( perl -wpe 'chomp; s/^/ --input /g' ./bams_inputs.list ) \
            --alleles ~{recallIndexedVcf.file} \
            --intervals ~{recallIndexedVcf.file} \
            --force-call-filtered-alleles \
            --force-active \
            --create-output-variant-index \
            --max-assembly-region-size 400 \
            --max-reads-per-alignment-start 200 \
            --output "~{outputVcfBasename}~{vcfSuffix}" \
            -contamination ~{default=0 contamination} \
            ~{ if useSpanningEventGenotyping then "" else " --disable-spanning-event-genotyping "} \
            ~{bamoutArg}

        touch ~{outputVcfBasename}.bamout.bam
    >>>
    
    output {
        File vcf = outputVcfBasename + vcfSuffix
        File vcfIdx = outputVcfBasename + vcfSuffix + ".tbi"
        IndexedFile vcfOut = {
          "file" : outputVcfBasename + vcfSuffix,
          "index" : outputVcfBasename + vcfSuffix + ".tbi"
        }
        IndexedFile idxVcf = {
          "file" : outputVcfBasename + vcfSuffix,
          "index" : outputVcfBasename + vcfSuffix + ".tbi"
        }
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}
#
#((ml GATK/4.2.4.1-Java-8-LTS && gatk --java-options "-Xmx4g" HaplotypeCaller -R 
#$(realpath 949052/Homo_sapiens_assembly38.fasta) --input $(realpath MySampleS2_aligned_reads.bam) \
# --input $(realpath MySampleS2_aligned_reads.bam) --input $(realpath  MySample ) \
# --alleles  variants_merged_scat0.vcf.gz --force-call-filtered-alleles \
# --force-active --create-output-variant-index --max-assembly-region-size 400 \
# --output test.vcf -L variants_merged_scat0.vcf.gz; ) 

task VariantEval {
    input {
        Array [File] inputVcfs
        Array[IndexedFile] inputIdxVcfs
        Array [String] inputTags
        File? KnownVariants #dbsnp
        Array [File] inputCompVcfs = inputVcfs
        Array [IndexedFile] inputCompIdxVcfs = inputIdxVcfs
        Array [String] inputCompTags = inputTags
        #File inputBai
        #File targetIntervalList
        String outputBasename
        Reference reference
        String gatkModule = "GATK"
        Int memoryGb = "4"
        Int javaMemoryGb = memoryGb - 1
        Int timeMinutes = 10 + ceil(size(inputVcfs, "G")) * 20
        Int? javaXmxMemoryMb = floor(memoryGb*0.9*1024)
        Int disk = ceil(size(inputVcfs, "M")*1.2)
    }
    #https://gatk.broadinstitute.org/hc/en-us/articles/30332019053723-VariantEval-BETA
    # gatk VariantEval \
    #-R reference.fasta \
    #-O output.eval.grp \
    #--eval set1:set1.vcf \
    #--eval set2:set2.vcf \
    #[--comp comp.vcf]
 
    command {
        set -e
        EVALSETS=""
        while IFS="" read -r LINE; do
            EVALSETS="$EVALSETS"" --eval:$LINE"
        done < <(paste ~{write_lines(select_all(inputTags))} ~{write_lines(select_all(inputVcfs))} -d \  )
        COMPSETS=""
        while IFS="" read -r LINE; do
            EVALSETS="$EVALSETS"" --comp:$LINE"
        done < <(paste ~{write_lines(select_all(inputCompTags))} ~{write_lines(select_all(inputCompVcfs))} -d \  )

        ml ~{gatkModule}
        gatk --java-options "-Xmx${javaXmxMemoryMb}m -Xms${javaXmxMemoryMb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            VariantEval \
            -R ~{reference.fasta} \
            $EVALSETS \
            $COMPSETS \
            -O "~{outputBasename}.grp" || touch "~{outputBasename}.grp"

        
    }
    
    output {
        File comparison = outputBasename + '.grp'
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}