version 1.0

task MultiQC {
    input {
        #should be files or dirs
        Array[File] files
        Array[File?] optionalFiles
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
        cat ~{write_lines(select_all(optionalFiles))} ~{write_lines(files)} | \
            (while read FILE; do 
                if [ -e "$FILE" ]; then
                    if [[ $FILE == *.zip ]]; then
                        mkdir -p "./unzip/$(basename $FILE .zip)"
                        ( cd "./unzip/$(basename $FILE .zip)" && unzip "$FILE" &>>./unzip.log)
                        find "./unzip/$(basename $FILE .zip)" -type f
                    else
                        echo "$FILE"
                    fi 
                fi
        done ) > ./filelist.txt
        module load ~{multiqcModule} && multiqc --force --file-list ./filelist.txt
    >>>

    output {
        File dir =  "multiqc_data"
        File html = "multiqc_report.html"
    }

    runtime {
        memory: select_first([memoryGb * 1024,1024])
        timeMinutes: timeMinutes
    }
}