version 1.0
import "../structs.wdl" as structs


task Callerise {
    input {
        #query sorted unaligned bam
        File inputVcf
        String outputBase
        String caller
        Int memoryGb = "1"
        #Int javaXmxMemoryMb = floor((memoryGb-0.5)*0.95*1024)
        String pipelineUtilModule = "pipeline-util"
        Int timeMinutes = 5 + ceil(size(select_all([inputVcf]), "G")) * 50
    }
    command {
        #also note the optional "-p" to only retain the tag prefixed data 
        
        module load ~{pipelineUtilModule} && \
        CalleriseVcf.pl -c ~{caller} -i ~{inputVcf} -o /dev/stdout | bgzip -c > "~{outputBase}.vcf.gz"
        tabix -p vcf "~{outputBase}.vcf.gz"
        
    }

    output {
        File vcf = outputBase+".vcf.gz"
        File vcfIdx = outputBase+".vcf.gz.tbi"
        
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
    }
}
task ReannotateVariants {
    input{
        IndexedFile combinedVariants
        Array [File] inputVcfsFiles
        Array [IndexedFile] inputVcfs
        String outputBasename
        Int memoryGb = "1"
        String pipelineUtilModule = "pipeline-util"
        String vcfSuffix = ".vcf.gz"
        Int timeMinutes = 1 + ceil(size(inputVcfsFiles, "G")) * 40
        #Possible values: {unsorted, queryname, coordinate, duplicate, unknown} #
        Int disk = 1 + ceil(size(inputVcfsFiles, "G")) * 1024
        
    }
    command <<<
        set -eo pipefail
        ml ~{pipelineUtilModule}
        if [ $(gzip -qdc  ~{combinedVariants.file} | grep -vc '^#' ) -eq 0 ]; then 
            if [[ "~{combinedVariants.file}" == *.vcf.gz ]]; then
                >&2 echo " ## "$(date)" ## Only header file detected with known extension copying input to output" 
                cp $(realpath ~{combinedVariants.file}) "~{outputBasename}~{vcfSuffix}"
                cp $(realpath ~{combinedVariants.index}) "~{outputBasename}~{vcfSuffix}.tbi"
            else
                >&2 echo " ## "$(date)" ## ERROR: Only header file detected with unknown extension. Exiting" 
                exit 1 
            fi
        else
            >&2 echo " ## "$(date)" ## Reannotating with older data." 
            perl $EBROOTPIPELINEMINUTIL/bin/RecoverSampleAnnotationsAfterCombineVariantsByPosWalk.pl \
                ~{outputBasename}.complex.vcf \
                ~{combinedVariants.file} \
                ~{sep=' ' inputVcfsFiles} \
                |bgzip -c >  ~{outputBasename}~{vcfSuffix}
            
            tabix -p vcf ~{outputBasename}~{vcfSuffix}
        fi
    >>>
    output {
        File vcf = outputBasename + vcfSuffix
        File vcfIdx = vcf + ".tbi"
        IndexedFile vcfOut = { 
          "file" : vcf,
          "index" : vcfIdx
        }
        IndexedFile idxVcf = { 
          "file" : vcf,
          "index" : vcfIdx
        }
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}
task AdFilter {
    input{
        IndexedFile inputVariantsToFilter
        Array [File] inputVcfsFiles
        Array [IndexedFile] inputVcfs
        String outputBasename
        Float freq = 0.01
        Int count = 4
        Float normalFreq = 0.2
        Int memoryGb = "1"
        String pipelineUtilModule = "pipeline-util"
        String vcfSuffix = ".vcf.gz"
        Int timeMinutes = 1 + ceil(size(inputVcfsFiles, "G")) * 40
        Int disk = 1 + ceil(size(inputVcfsFiles, "G")) * 1024
        
    }
    #this filters for allele depth control(s) + 4 and/or frequency control(s) + 0.05. 
    #Maybe also omit entirely when seen in controls above a certain freq like 0.10 pct.
    command <<<
        set -eo pipefail
        ml ~{pipelineUtilModule}
        if [ -e ~{inputVariantsToFilter.file} ] && \
            [  "$(bgzip -cd ~{inputVariantsToFilter.file} | grep -vc '^#')" != 0 ]; then
            perl $EBROOTPIPELINEMINUTIL/bin/AdFilter.pl \
             ~{"-f " + freq} \
             ~{"-n " + normalFreq} \
             ~{"-c " + count} \
            "~{inputVariantsToFilter.file}" \
            "~{sep='" "' inputVcfsFiles}" | \
             perl -wane 'print if(m/^#/ || (
                m/\tPASS\t/ && (( defined($F[9]) && $F[9] ne "." ) || defined($F[11]) ) ));' | \
            bgzip -c >  "~{outputBasename}~{vcfSuffix}"
        else
            if [[ "~{inputVariantsToFilter.file}" == *.gz ]]; then
                cp "~{inputVariantsToFilter.file}"  "~{outputBasename}~{vcfSuffix}"
            else
               bgzip -c "~{inputVariantsToFilter.file}" > "~{outputBasename}~{vcfSuffix}"
            fi
        fi
        
        tabix -p vcf "~{outputBasename}~{vcfSuffix}"

    >>>
    output {
        File vcf = outputBasename + vcfSuffix
        File vcfIdx = vcf + ".tbi"
        IndexedFile vcfOut = { 
          "file" : vcf,
          "index" : vcfIdx
        }
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}