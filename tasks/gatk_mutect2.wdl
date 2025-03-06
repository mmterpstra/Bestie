version 1.0

import "../structs.wdl" as structs

task MuTect2JointCalling {
    input {
        #Array[IndexedFile] inputIndexedBams
        Array[File] inputBams
        Array[File] inputBamIndexes
        Array[String] normals
        File targetIntervalList
        String outputVcfBasename
        String vcfSuffix = ".vcf.gz"
        Reference reference
        IndexedFile germlineAfVcf #aka gnomadOnlyAfVcf
        IndexedFile panelOfNormalsVcf
        #f1r2TarGz should probs be true default except on pon creation (not sure if this gives an performance penalty though)
        Boolean createF1r2TarGz = true
        #maxMnpDistance usually set to 0 for genomicsDBImport on pon creation default 1 to merge adjecent mnps
        Int? maxMnpDistance
        #maxReadsPerAlignmentStart controls downsampling default 50 0=disable
        Int? maxReadsPerAlignmentStart
        #default settings
        String gatkModule = "GATK"
        Int memoryGb = "8"
        Int javaMemoryGb = memoryGb - 1
        
        Int? javaXmxMemoryMb = floor(memoryGb*0.9*1024)
        Int disk = ceil(size(inputBams,
            "M")*1.2 + 5000)
        # it bugs out on ? structs
        #select_all([germlineAfVcf.file,germlineAfVcf.index]),
        #select_all([panelOfNormalsVcf.file,panelOfNormalsVcf.index])

        Int timeMinutes = ceil(10 + size(inputBams,
            "G")*1.2 ) * 50
    }
    
    #IndexedFile controlBam = select_first([inputControlBam,inputBam])
    #String normalSpec = if (artifactDetection) then " --artifact_detection_mode " else " -I:normal " + inputControlBam.file
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
        set -e
        set -o pipefail 

        ml ~{gatkModule}
        echo -n  "">"./bams_inputs.list"
        cat ~{write_lines(select_all(inputBams))} \
            ~{write_lines(inputBamIndexes)} \
            ~{write_lines(select_all([germlineAfVcf.file,germlineAfVcf.index]))}\
            ~{write_lines(select_all([panelOfNormalsVcf.file,panelOfNormalsVcf.index]))} | \
            (while read FILE; do 
                if [ ! -e "$TMPDIR/""$(basename "$FILE")" ]; then
                    cp "$FILE" "$TMPDIR/"
                    if [[ $FILE =~ \.bam$ ]]; then
                        echo " -I ""$TMPDIR/""$(basename "$FILE")"" ">>"./bams_inputs.list" 
                    fi
                else
                    echo "Duplicate file basename spotted $FILE" && exit 1 
                fi
            done )

        gatk --java-options "-Xmx~{javaXmxMemoryMb}m -Xms~{javaXmxMemoryMb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            Mutect2 \
            -R ~{reference.fasta} \
            $(cat "./bams_inputs.list" | python3 -c "import sys; print(' '.join(str(x.rstrip()) for x in sys.stdin))") \
            -normal ~{sep=' -normal ' normals} \
            ~{"--germline-resource $TMPDIR/$(basename " + germlineAfVcf.file + ") " }  \
            ~{"--panel-of-normals $TMPDIR/$(basename " + panelOfNormalsVcf.file + ") " } \
            ~{if createF1r2TarGz then "--f1r2-tar-gz " + outputVcfBasename + ".f1r2.tar.gz" else "" } \
            ~{"--max-mnp-distance " + maxMnpDistance } \
            ~{"--max-reads-per-alignment-start " + maxReadsPerAlignmentStart} \
            -L ~{targetIntervalList} \
            -O ~{outputVcfBasename}.vcf.gz

        ~{if ! createF1r2TarGz then "touch " + outputVcfBasename + ".f1r2.tar.gz.skipped" else "#didnt skip .f1r2.tar.gz" }

    }
    
    output {
        File vcf = outputVcfBasename + '.vcf.gz'
        File vcfIdx = outputVcfBasename +  '.vcf.gz.tbi'
        File stats = outputVcfBasename + '.vcf.gz.stats'
        File f1r2TarGz = if createF1r2TarGz then outputVcfBasename + '.f1r2.tar.gz' else outputVcfBasename + ".f1r2.tar.gz.skipped"
        IndexedFile vcfOut = {
          "file" : outputVcfBasename + '.vcf.gz',
          "index" : outputVcfBasename + '.vcf.gz.tbi'
        }
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}
task LearnReadOrientationModel {
    input {
        Array[File] f1r2TarGz
        String outputArtifactPriorTableBasename
        String gatkModule
        Int memoryGb = "8"
        Int javaMemoryGb = memoryGb - 1
        Int? javaXmxMemoryMb = floor(memoryGb*0.9*1024)
        Int disk = ceil(size(f1r2TarGz, "M")*1.2)
        Int timeMinutes = 1 + ceil(size(f1r2TarGz, "G")) * 30
    }
    command {
        set -e
        set -o pipefail 
        
        ml ~{gatkModule}
        gatk --java-options "-Xmx${javaXmxMemoryMb}m -Xms${javaXmxMemoryMb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            LearnReadOrientationModel \
            --input ~{sep=' --input ' f1r2TarGz} \
            --output "~{outputArtifactPriorTableBasename}"".tar.gz"
    }
    output {
        File artifactpriortable = outputArtifactPriorTableBasename + '.tar.gz'
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}
task MergeMutectStats {
    input {
        Array[File] inputMutectStats
        String outputMergedStats = "merged.stats"
        String gatkModule
        Int memoryGb = "2"
        Int javaMemoryGb = memoryGb - 1
        Int? javaXmxMemoryMb = floor(memoryGb*0.9*1024)
        Int disk = ceil(size(inputMutectStats, "M")*1024*1.2)
        Int timeMinutes = 1 + ceil(size(inputMutectStats, "G")) * 30
    }
    command {
        set -e

        ml ~{gatkModule}

        gatk --java-options "-Xmx${javaXmxMemoryMb}m -Xms${javaXmxMemoryMb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            MergeMutectStats \
            -stats ~{sep=' -stats ' inputMutectStats} \
            --output "~{outputMergedStats}"
    }
    output {
        File stats = outputMergedStats
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}
task FilterMutect {
    input {
        IndexedFile inputSomaticVcf
        File? contaminationTable
        File? segmentsTsv
        File artifactPriorsTarGz
        File? stats 
        String vcfSuffix = ".vcf.gz"
        File targetIntervalList
        String outputVcfBasename
        Reference reference
        String gatkModule = "GATK"
        Int memoryGb = "4"
        Int javaMemoryGb = memoryGb - 1
        Float? contamination = 0
        Int timeMinutes = 1 + ceil(size(inputSomaticVcf.file, "G")) * 120
        Int? javaXmxMemoryMb = floor(memoryGb*0.9*1024)
        Int disk = ceil(size(inputSomaticVcf.file, "M")*1.2)
    }

    command {
        set -e

        ml ~{gatkModule}
        gatk --java-options "-Xmx${javaXmxMemoryMb}m -Xms${javaXmxMemoryMb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
                FilterMutectCalls \
                -R ~{reference.fasta} \
                -V ~{inputSomaticVcf.file} \
                ~{"--contamination-table " + contaminationTable} \
                ~{"--tumor-segmentation " + segmentsTsv} \
                ~{"-stats " + stats} \
                --ob-priors ~{artifactPriorsTarGz} \
                --output ~{outputVcfBasename}~{vcfSuffix}
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