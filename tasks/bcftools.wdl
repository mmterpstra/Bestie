version 1.0

import "../structs.wdl" as structs


task Index {
    input {
        File inputVcf
        Int memoryGb = "1"
        String bcftoolsModule = "BCFtools/1.21-GCC-12.2.0"
        Int timeMinutes = 1 + ceil(size(inputVcf, "G")) * 120
        #Possible values: {unsorted, queryname, coordinate, duplicate, unknown} #
        Int disk = 1 + ceil(size(inputVcf, "G")) * 1024

    }
    command {
        set -e
        #cp ~{inputVcf} ./
        module load ~{bcftoolsModule} && \
        bcftools index -t ~{inputVcf}
        #this should make it 100% findable by cromwell
        cp "~{inputVcf}" "~{inputVcf}.tbi" ./
    }

    output {
        File outputVcf = inputVcf
        File outputVcfIdx = inputVcf + ".tbi"
        IndexedFile vcfOut = {
          "file" : inputVcf,
          "index" : inputVcf + ".tbi"
        }
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}
