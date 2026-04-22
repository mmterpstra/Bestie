version 1.0
import "../structs.wdl"

task MultiQC {
    input {
        #should be files or dirs
        Array[File] files
        Array[File?] optionalFiles
        Array[String] configSettings = ["config settings"]
        File moduleLog
        SampleConfig sampleConfig
        String prefix = "multiqc"
        #Scales badly might be a multiple of files
        Int? memoryGb = 8
	    String multiqcModule
        Int timeMinutes = 120
    }

    command <<<
        set -e


        (
            #yaml creation
            echo "title: 'Multiqc report of FastqToBam workflow'"
            echo "subtitle: 'Rundate $(date) workflow version $VERSION' "
            echo "intro_text: 'The conversion of the fastq files to variant call ready bam files with tools, versions and settings shown below' "
            echo 'report_header_info:' 
            echo "- 'Used modules' : ''" 
            cat ~{moduleLog} | python3 -c "import sys; [ print('- \'\': \''+line.rstrip()+'\'' ) for count,line in enumerate(sys.stdin)]"  
            echo "- '' : ''" 
            echo "- 'Used settings' : ''" 
            cat ~{write_lines(select_all(configSettings))}| \
                grep -v 'Testing\|----\|Where\|Default Module' |\
                python3 -c "import sys; [ print('- \'\': \''+line.rstrip()+'\'' ) for count,line in enumerate(sys.stdin)]"  
            echo "- '' : ''" 
            echo "- Sample settings : \|"
            cat <<- SAMPLEJSONHEADER
            This config was used for the sample specific settings.

            ```json
            SAMPLEJSONHEADER
            cat ~{write_json(sampleConfig)}
            echo '```'
            echo '``'

            echo "- '': 'For details about the used tools reference the readme https://github.com/mmterpstra/Bestie/blob/develop/README.md'"

            echo "- Readme : ''"
            echo "- '': 'For details about the used tools reference the readme https://github.com/mmterpstra/Bestie/blob/develop/README.md'"

        ) >  "multiqc.yaml"
        #cat ~{write_lines(files)} > filelist.txt
        #test the optional files , unzip if needed and add to filelist 
        #This is due to multiqc not handling zip archives
        #https://github.com/ewels/MultiQC/blob/f4878540c0b1e77fa419b589ec33f3a5039d17da/multiqc/modules/picard/HsMetrics.py#L92-L96 it looks at the bam "basename" cleaned up to determine witch sample it is derived from.
        #this skips localisation of optional files so wip
        DIRNO=1
        cat ~{write_lines(select_all(optionalFiles))} ~{write_lines(files)} | \
            (while read FILE; do 
                if [ -e "$FILE" ]; then
                    if [[ $FILE == *.zip ]]; then
                        mkdir -p "./unzip/${DIRNO}_$(basename $FILE .zip)"
                        ( cd "./unzip/${DIRNO}_$(basename $FILE .zip)" && unzip "$FILE" &>>./unzip.log) || \
                            (>&2 echo "## ERROR ## 'cd $PWD/unzip/${DIRNO}_$(basename $FILE .zip)' or 'unzip $FILE' failed log in './unzip.log'." && exit 1)
                        find "./unzip/${DIRNO}_$(basename $FILE .zip)" -type f
                    else
                        echo "$FILE"
                    fi 
                fi
                ((DIRNO++))
        done ) > ./filelist.txt
        module load ~{multiqcModule} && multiqc --force --filename ~{prefix} --file-list ./filelist.txt -c multiqc.yaml
    >>>

    output {
        File dir =  prefix+ "_data"
        File html = prefix + ".html"
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
    }
}