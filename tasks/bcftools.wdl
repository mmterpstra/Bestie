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
        set -e -o pipefail
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
        Int disk = 1 + ceil(size(inputVcf, "G")) * 1024
    }
    command {
        set -e -o pipefail
        #cp ~{inputVcf} ./
        module load ~{bcftoolsModule} && \
        bcftools norm \
            --atomize \
            --atom-overlaps . \
            --fasta-ref ~{reference.fasta} \
            --old-rec-tag OLD_RECORD \
            --write-index=tbi \
            -m -any \
            --sort lex \
            --output ~{outputBasename}".vcf.gz" \
            ~{inputVcf} 

    }

    output {
        File vcf = outputBasename + ".vcf.gz"
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
        perl -wpe 's/\t\.:\.(:\.)+/\t./ if m/\t\.:\.(:\.|:[ATCG]+)+[\t\n]/' | \
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
        IndexedFile idxVcf = {
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
    #bcftools isec with the twist that it spits out a really
    # basic but bcftools compatible vcf file for more tooling options 
    input {
        Array [File] inputVcfs
        Array [IndexedFile] inputIdxVcfs
        String outputBasename
        Int memoryGb = "1"
        String bcftoolsModule = "BCFtools/1.21-GCC-12.2.0"
        Boolean applyFilters = false
        String filters = 'PASS,.'
        Int timeMinutes = 1 + ceil(size(inputVcfs, "G")) * 120
        #Possible values: {unsorted, queryname, coordinate, duplicate, unknown} #
        Int disk = 1 + ceil(size(inputVcfs, "G")) * 1024
        Int minIntersect = 1

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
        set -e -o pipefail
        module load ~{bcftoolsModule} 
        NFILES=$(ls ~{sep=' ' inputVcfs} | wc -l )
        if [ $NFILES -ge 2 ] && \
            [ $(gzip -qdc ~{sep=' ' inputVcfs} | head -n 50000 | grep -cv '^#') -gt 0 ]; then


            (
                bgzip -dc  ~{inputVcfs[0]} | grep '##fileformat\|##contig' ;
                bcftools isec \
                    -c all --nfiles +~{minIntersect} \
                    --output /dev/stdout \
                    ~{sep=' ' inputVcfs} \
                    ~{if applyFilters then "--apply-filters "+ filters else "" } | \
                    python3 -c "import sys; print( \
                    '##INFO=<ID=ISEC,Number=1,Type=String,Description=\"bcftools isec results as a 0/1 presence table from files ~{sep=' ' inputVcfs}\">\n' + \
                    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'); \
                    [print('\t'.join([fields[0],fields[1],'.',fields[2],fields[3],'.','.','ISEC='+fields[4]])) for line in sys.stdin if (fields:=line.rstrip('\n').split('\t'))]"
            )| bgzip -c > ~{outputBasename}.vcf.gz
            tabix -p vcf ~{outputBasename}.vcf.gz

        else
            module load ~{bcftoolsModule} 

            cat ~{sep=' ' inputVcfs} >  ~{outputBasename}.vcf.gz
            tabix -p vcf ~{outputBasename}.vcf.gz
        fi
    >>>
    

    output {
        File vcf = outputBasename + ".vcf.gz"
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
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}

task Concat {
    input {
        Array [File] inputVcfs
        #the following is not used explictly but implicitly by using `find ../inputs/*`
        Array [IndexedFile] inputIndexedVcfs

        String outputBasename
        Int memoryGb = "1"
        String bcftoolsModule = "BCFtools/1.21-GCC-12.2.0"
        Int timeMinutes = 1 + ceil(size(inputVcfs, "G")) * 120
        #Possible values: {unsorted, queryname, coordinate, duplicate, unknown} #
        Int disk = 1 + ceil(size(inputVcfs, "G")) * 1024

    }
    command {
        set -e -o pipefail

        >&2 echo " ## "$(date)" ## Localising files"
        #finds both files and links due to how the linking system can work. 
        /usr/bin/find ../inputs/* -type l -o -type f | \
            (while read FILE; do 
                if [ ! -e "$TMPDIR/""$(basename "$FILE")" ]; then
                    cp "$FILE" "$TMPDIR/"
                    if [[ $FILE =~ \.vcf.gz$ ]]; then
                        echo "$TMPDIR/""$(basename "$FILE")">>"./vcfs_inputs.list" 
                    fi
                else
                    echo "Duplicate file basename spotted $FILE" && exit 1 
                fi
            done )
        >&2 echo " ## "$(date)" ## Loading modules" 
        module load ~{bcftoolsModule}
        
        >&2 echo " ## "$(date)" ## Running bcftools " 

        #the perl  removes .:.:.:.:.:.:.:.:. bs in the sample descriptions
        #ls  "~{sep="\" \"" inputVcfs}"
        #this catches the single scatter interval segfault of bcftools 
        if [ "$(wc -l "./vcfs_inputs.list" )" -eq 2 ]; then
            (cat "./vcfs_inputs.list")  | \
            (
                while read FILE; do 
                    if [[ $FILE =~ \.vcf.gz$ ]]; then
                        cp "$FILE" ~{outputBasename}".vcf.gz"
                        cp  "$FILE"".tbi"  ~{outputBasename}".vcf.gz.tbi"
                    fi
                done
            )
        else 
        
            bcftools concat \
            --allow-overlaps \
            --rm-dups exact \
            $(cat ./vcfs_inputs.list) | \
            perl -wpe 's/\t\.:\.(:\.)+/\t./ if m/\t\.:\.(:\.)+[\t\n]/' | \
            bgzip -c > ~{outputBasename}".vcf.gz"
            tabix -p vcf  ~{outputBasename}".vcf.gz"

        fi
    }

    output {
        File vcf = outputBasename + ".vcf.gz"
        File vcfIdx = vcf + ".tbi"
        IndexedFile idxVcf = {
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
#task GtCheck {
#
#}
task Stats {
    input {
        Array [File] inputVcfs
        Array [IndexedFile] inputIndexedVcfs
        String outputBasename
        Int memoryGb = "1"
        String bcftoolsModule = "BCFtools/1.21-GCC-12.2.0"
        Int timeMinutes = 1 + ceil(size(inputVcfs, "G")) * 120
        #Possible values: {unsorted, queryname, coordinate, duplicate, unknown} #
        Int disk = 1 + ceil(size(inputVcfs, "G")) * 1024

    }
    command {
        set -e 
        set -o pipefail

        >&2 echo " ## "$(date)" ## Localising files"
        #finds both files and links due to how the linking system can work. 
        /usr/bin/find ../inputs/* -type l -o -type f | \
            (while read FILE; do 
                if [ ! -e "$TMPDIR/""$(basename "$FILE")" ]; then
                    cp "$FILE" "$TMPDIR/"
                    if [[ $FILE =~ \.vcf.gz$ ]]; then
                        echo "$TMPDIR/""$(basename "$FILE")">>"./vcfs_inputs.list" 
                    fi
                else
                    echo "Duplicate file basename spotted $FILE" && exit 1 
                fi
            done )
        >&2 echo " ## "$(date)" ## Loading modules" 
        module load ~{bcftoolsModule}
        
        >&2 echo " ## "$(date)" ## Running bcftools " 

        #this catches the single scatter interval segfault of bcftools 
        if [ "$(wc -l < "./vcfs_inputs.list" )" -eq 1 ]; then
            while IFS= read -r FILE; do
                if [[ $FILE == *.vcf.gz ]]; then
                    bcftools stats "$FILE" > "~{outputBasename}"".bcftools_stats"
                else 
                    >&2 echo "## ""$(date)"" ## Warning ## $FILE is not recognised as .vcf.gz file. "
                fi
            done < "./vcfs_inputs.list"
        else 
            #this fifo here is to improve multiqc output 
            # Since it tends to only take a single file with the same output name 

            mkfifo  ~{outputBasename}.tmp.vcf

            bcftools concat \
            --allow-overlaps \
            --rm-dups exact \
            $(cat ./vcfs_inputs.list) >  ~{outputBasename}.tmp.vcf &
            
            bcftools stats --af-bins 0.005,0.01,0.02,0.05,0.10,0.20,0.50,0.80,1 ~{outputBasename}.tmp.vcf  > ~{outputBasename}".bcftools_stats"

        fi
    }

    output {
        File stats = outputBasename + ".bcftools_stats"
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}

task Convert {
    input {
        IndexedFile inputIndexedVcf
        Array [String] sampleNames
        String outputBasename
        Int memoryGb = "1"
        String bcftoolsModule = "BCFtools/1.21-GCC-12.2.0"
        Int timeMinutes = 5 + ceil(size(inputIndexedVcf.file, "G")) * 20
        #Possible values: {unsorted, queryname, coordinate, duplicate, unknown} #
        Int disk = 1 + ceil(size(inputIndexedVcf.file, "G")) * 1024

    }
    command {
        set -e -o pipefail

        >&2 echo " ## "$(date)" ## Localising files"
        
        module load ~{bcftoolsModule}
        
        >&2 echo " ## "$(date)" ## Running bcftools " 

       #the --samples doesnt seem to work in convert idk why
        bcftools convert \
        --samples '~{sep=',' sampleNames}' \
        ~{inputIndexedVcf.file} | \
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