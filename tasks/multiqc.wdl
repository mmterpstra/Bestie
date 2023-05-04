version 1.0

task MultiQC {
    input {
        #should be files or dirs
        Array[File] files
        Array[String] optionalFiles
        #Scales badly
        Int? memoryGb = "4"
	    String multiqcModule
        Int timeMinutes = 50
    }

    command <<<
        set -e
        cat ~{write_lines(files)} > filelist.txt
        (while read FILE; do 
            if [ -e $FILE ]; then
                echo "$FILE"
            fi
        done <  ~{write_lines(optionalFiles)}) >> ./filelist.txt
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