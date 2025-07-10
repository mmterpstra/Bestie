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
task Norm {
    input {
        File inputVcf
        String outputBasename
        Reference reference
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
        bcftools norm \
            --atomize \
            --atom-overlaps . \
            --fasta-ref ~{reference.fasta} \
            --old-rec-tag OLD_RECORD \
            --write-index=tbi \
            -m -any \
            --output ~{outputBasename}".vcf.gz" \
            ~{inputVcf} \

    }

    output {
        File vcf = outputBasename + ".vcf.gz"
        File vcfIdx = vcf + ".tbi"
        IndexedFile vcfOut = {
          "file" : vcf,
          "index" : vcfIdx
        }
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}

task ViewSamples {
    input {
        File inputVcf
        Array [String] samples
        String outputBasename
        Int memoryGb = "1"
        String bcftoolsModule = "BCFtools/1.21-GCC-12.2.0"
        Int timeMinutes = 1 + ceil(size(inputVcf, "G")) * 120
        #Possible values: {unsorted, queryname, coordinate, duplicate, unknown} #
        Int disk = 1 + ceil(size(inputVcf, "G")) * 1024

    }
    command {
        set -eo pipefail
        #cp ~{inputVcf} ./
        #the perl  removes .:.:.:.:.:.:.:.:. bs in the sample descriptions
        module load ~{bcftoolsModule} && \
        bcftools view \
            --samples ~{sep=',' samples} \
            ~{inputVcf} | \
        perl -wpe 's/\t\.:\.(:\.)+/\t./ if m/\t\.:\.(:\.)+[\t\n]/' | \
        bgzip -c > ~{outputBasename}".vcf.gz"
        tabix -p vcf  ~{outputBasename}".vcf.gz"
    }

    output {
        File vcf = outputBasename + ".vcf.gz"
        File vcfIdx = vcf + ".tbi"
        IndexedFile vcfOut = {
          "file" : vcf,
          "index" : vcfIdx
        }
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}

task Isec {
    input {
        Array [File] inputVcfsFiles
        Array [IndexedFile] inputVcfs
        String outputBasename
        Int memoryGb = "1"
        String bcftoolsModule = "BCFtools/1.21-GCC-12.2.0"
        Int timeMinutes = 1 + ceil(size(inputVcfsFiles, "G")) * 120
        #Possible values: {unsorted, queryname, coordinate, duplicate, unknown} #
        Int disk = 1 + ceil(size(inputVcfsFiles, "G")) * 1024

    }
    #This outputs an isec specific table
    #chr8:127736594-127740958        2501    A       G       11
    #chr8:127736594-127740958        3401    G       A       11
    #chr8:127736594-127740958        3701    T       C       11
    #chr8:127736594-127740958        3821    G       A       11
    #chr8:127736594-127740958        4145    G       A       10
    #chr8:127736594-127740958        4346    A       AT      10
    #chr8:127736594-127740958        4363    A       T       10

    #reformatting the table:
    #CHROM  POS REF ALT BINARYPRESENCE_TABLE  
    #to (omitting the BINARYPRESENCE_TABLE column for now)
    ##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
    command <<<
        set -e
        module load ~{bcftoolsModule} && \
        bcftools isec \
            -c all --nfiles +1 \
            --output /dev/stdout \
            ~{sep=' ' inputVcfsFiles} | \
            python3 -c "import sys; print('##fileformat=VCFv4.2\n'+ \
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'); \
            [print('\t'.join([fields[0],fields[1],'.',fields[2],fields[3],'.','.','.'])) for line in sys.stdin if (fields:=line.rstrip('\n').split('\t'))]"| \
            bgzip -c > ~{outputBasename}.vcf.gz
        tabix -p vcf ~{outputBasename}.vcf.gz
    >>>
    

    output {
        File vcf = outputBasename + ".vcf.gz"
        File vcfIdx = vcf + ".tbi"
        IndexedFile vcfOut = {
          "file" : vcf,
          "index" : vcfIdx
        }
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}

task Concat {
    input {
        Array [File] inputVcfs
        Array [IndexedFile] inputVcfsIndexed
        String outputBasename
        Int memoryGb = "1"
        String bcftoolsModule = "BCFtools/1.21-GCC-12.2.0"
        Int timeMinutes = 1 + ceil(size(inputVcf, "G")) * 120
        #Possible values: {unsorted, queryname, coordinate, duplicate, unknown} #
        Int disk = 1 + ceil(size(inputVcf, "G")) * 1024

    }
    command {
        set -eo pipefail
        #the perl  removes .:.:.:.:.:.:.:.:. bs in the sample descriptions
        module load ~{bcftoolsModule} && \
        bcftools concat \
         --allow-overlaps \
         --rm-dups exact \
         "~{sep="\" \\\n\"" inputVcfs}" | \
        perl -wpe 's/\t\.:\.(:\.)+/\t./ if m/\t\.:\.(:\.)+[\t\n]/' | \
        bgzip -c > ~{outputBasename}".vcf.gz"
        tabix -p vcf  ~{outputBasename}".vcf.gz"
    }

    output {
        File vcf = outputBasename + ".vcf.gz"
        File vcfIdx = vcf + ".tbi"
        IndexedFile vcfOut = {
          "file" : vcf,
          "index" : vcfIdx
        }
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}