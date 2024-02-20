version 1.0

task MultiQC {
    input {
        #should be files or dirs
        Array[File] files
        Array[File?] optionalFiles
        String prefix = "multiqc"
        #Scales badly might be a multiple of files
        Int? memoryGb = 8
	    String multiqcModule
        Int timeMinutes = 120
    }

    command <<<
        set -e
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
        module load ~{multiqcModule} && multiqc --force --filename ~{prefix} --file-list ./filelist.txt
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