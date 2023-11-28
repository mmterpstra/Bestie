version 1.0
import "../structs.wdl" as structs


task CollapseFastq {
    input {
        Array[File] reads
        #String dockerImage = "quay.io/biocontainers/cutadapt:4.2--py310h1425a21_0"
        String outputPrefix
        String? threeLetterName = "AAA"
        Int memoryGb = 8
        Boolean nextflex = false
        Int timeMinutes = 10 + ceil(size(reads, "GiB") * 40 )
        String? dockerImage
        #String dockerImage = "quay.io/biocontainers/coreutils" #idk needs python3 and coreutils
    }
    #this onliner trims the first and last four bases from the sequences
    String nextFlexCmd = if (nextflex) then " python3 -c \"import sys;[sys.stdout.write(line.rstrip()[4:-4]+'\\n') for count,line in enumerate(sys.stdin)]\" |  " else ""

    #Python oneliner1 only prints the sequence string from the fastq.
    #Python oneliner2 assumes a stream of sequences delimited by newline prints the stuff if not newline 
    command <<<
        set -e -o pipefail
        mkdir -p "$(dirname ~{outputPrefix})"
        zcat ~{sep=" " reads} | \
        python3 -c "import sys; [ print(line) if count % 4 == 1 else None for count,line in enumerate(sys.stdin)]" | \
        ~{nextFlexCmd} \
        sort --temporary-directory="$(dirname ~{outputPrefix})" --buffer-size=4G | \
        uniq -c | \
        python3 -c "import sys;[None if(line.split(' ')[-1] == '\n') else sys.stdout.write('>'+sys.argv[1]+'_'+str(count)+'_x'+'\n'.join(line.split(' ')[-2:])) for count,line in enumerate(sys.stdin)]" ~{threeLetterName} > \
         "~{outputPrefix}"".md.fa"
    >>>
    output {
        File outputCollapsedFasta = outputPrefix + ".md.fa"
    }
    runtime {
        memory: memoryGb*1024
        timeMinutes: timeMinutes
        #docker: dockerImage
    }

}

#
#Borrowed stuff below here

#https://github.com/biowdl/tasks/blob/c92755e510723da731ba92637c41e58c8490b5fc/common.wdl#L66
# Copyright (c) 2017 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

task AppendToStringArray {
    input {
        Array[String] array
        String string
        
        Int memory = 256
    }

    command {
        echo "~{sep='\n' array}
        ~{string}"
    }

    output {
        Array[String] outArray = read_lines(stdout())
    }

    runtime {
        memory: memory
    }
}

task CreateLink {
    # Making this of type File will create a link to the copy of the file in
    # the execution folder, instead of the actual file.
    # This cannot be propperly call-cached or used within a container.
    input {
        String inputFile
        String outputPath

        Int memory = 256
    }

    command {
        echo $PWD
        ln -sf "~{inputFile}"  "~{outputPath}"
    }

    output {
        File link = outputPath
    }

    runtime {
        memory: memory
        timeMinutes: 10
    }
}

task CatSampleDescriptorJson {
    # This creates a new sampledescriptor from input sample descriptor later there should be ways to add in files/annotations
    input {
        SampleDescriptor sample
        Int memory = 256
    }
    command {
        cat ${write_json(sample)}
    }

    output {
        SampleDescriptor sampleUpdated = read_json(stdout())
    }

    runtime {
        memory: memory
        timeMinutes: 10
    }
}
task AddAlignedReadsToSampleDescriptor {
    # This creates a new sampledescriptor from input sample descriptor later there should be ways to add in files/annotations
    input {
        SampleDescriptor sample
        IndexedFile bam
        Int memory = 256
    }
    command <<<
        cat ~{write_json(sample)} | perl -wpe 'BEGIN{our $bam = shift(@ARGV);our $bamidx = shift(@ARGV);};s!"alignedReads":null!"alignedReads":{"file":"$bam","index":"$bamidx"}!g' "~{bam.file}" "~{bam.index}"    
    >>>

    output {
        SampleDescriptor sampleUpdated = read_json(stdout())
    }

    runtime {
        memory: memory
        timeMinutes: 10
    }
}

task ConcatenateTextFiles {
    input {
        Array[File] fileList
        String combinedFilePath
        Int memory = 256
    }

    # When input and output is both compressed decompression is not needed.
    #String cmdPrefix = if (unzip && !zip) then "zcat " else "cat "
    #String cmdSuffix = if (!unzip && zip) then " | gzip -c " else ""

    command {
        set -e -o pipefail
        mkdir -p "$(dirname ~{combinedFilePath})"
        cat ~{sep=" " fileList} > ~{combinedFilePath}
    }

    output {
        File combinedFile = combinedFilePath
    }

    runtime {
        memory: memory
    }
}
task ZipFiles {
    input {
        Array[File] fileList
        Array[String] optionalFileList = ["test.txt"]
        String outputPrefix
        Int memory = 512
        Boolean flattenArchive = true
    }
    
    #WIP: does not localise optional files...

    command {
        set -e -o pipefail
        cat ~{write_lines(fileList)} ~{write_lines(select_all(optionalFileList))} | 
        ( while read FILE; do 
            if [ -e $FILE ]; then
                echo "$FILE"
            fi
        done ) > ./filelist.txt

        pwd
        #try find cromwell execution dir and execute from there
        EXECUTION_DIR=$(head -n 1 ./filelist.txt | \
        python3 -c "import sys;[sys.stdout.write( line[:(line.find('/cromwell-executions/')+21)]+'\n') for line in sys.stdin]" )
        OUTDIR=$PWD

        #zip command
        cat ./filelist.txt | \
            python3 -c "import sys;[sys.stdout.write( line[(line.find('/cromwell-executions/')+21):]+'\n') for line in sys.stdin]" | \
            (cd $EXECUTION_DIR && \
                zip ~{if flattenArchive then "-j " else ""}  $OUTDIR/~{outputPrefix}.zip -@
            )
    }

    output {
        File zip = outputPrefix + ".zip"
    }

    runtime {
        memory: memory
    }
}