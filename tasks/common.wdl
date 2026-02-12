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
        ~{string}" > array_string.list
    }

    output {
        Array[String] outArray = read_lines("./array_string.list")
    }

    runtime {
        memory: memory
        timeMinutes: 20
    }
}

task UniqueArray {
    input {
        Array[String] array        
        Int memory = 256
    }

    command {
        echo -en "~{sep='\n' array}"| sort -u > ./unique.list
    }

    output {
        Array[String] outArray = read_lines("./unique.list")
    }

    runtime {
        memory: memory
        timeMinutes: 20
    }
}

task InverseSelection {
    input {
        Array[String] array
        Array[String] selection        
        Int memory = 256
    }

    command {
        #all data is 'com'pared and only lines unique to the first input file are emitted
        comm -23 <(sort -u ~{write_lines(array)}) <(sort -u ~{write_lines(selection)}) \
             > ./result.list
    }

    output {
        Array[String] result = read_lines("./result.list")
    }

    runtime {
        memory: memory
        timeMinutes:5
    }
}

task CreateIndexedLink {
    # Making this of type File will create a link to the copy of the file in
    # the execution folder, instead of the actual file.
    # This cannot be propperly call-cached or used within a container.
    input {
        File inputFile
        File indexFile
        Boolean hardLink = false
        Int memory = 256
    }

    command {
        echo $PWD 
        mkdir -p index/
        mkdir -p file/
        ln -t "./file/" -~{if hardLink then "" else "s" }f "$(realpath "~{inputFile}")"  
        ln -t "./index/" -~{if hardLink then "" else "s" }f "$(realpath "~{indexFile}")"  
    }

    output {
        File file = select_first(glob( "./file/*"))
        File index = select_first(glob( "./index/*"))
        IndexedFile out = {"file":file,"index":index}
    }

    runtime {
        memory: memory
        timeMinutes: 5
    }
}

task CreateLink {
    # Making this of type File will create a link to the copy of the file in
    # the execution folder, instead of the actual file.
    # This cannot be propperly call-cached or used within a container.
    input {
        String inputFile
        String outputPath
        Boolean hardLink = false
        Int memory = 256
    }

    command {
        echo $PWD
        ln -~{if hardLink then "" else "s" }f "~{inputFile}"  "~{outputPath}"
    }

    output {
        File link = outputPath
    }

    runtime {
        memory: memory
        timeMinutes: 5
    }
}
task CatSampleDescriptorJson {
    # This creates a new sampledescriptor from input sample descriptor later there should be ways to add in files/annotations
    input {
        SampleDescriptor sample
        Int memory = 256
    }
    command {
        cat ${write_json(sample)} > ~{sample.name}.json
    }

    output {
        SampleDescriptor sampleUpdated = read_json(sample.name + '.json')
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
        # this either adds alignedreads or replaces them.
        cat ~{write_json(sample)} | \
            perl -wpe 'BEGIN{our $bam = shift(@ARGV);our $bamidx = shift(@ARGV);};s!"alignedReads":null|"alignedReads"\:\{.*?\}!"alignedReads":{"file":"$bam","index":"$bamidx"}!g' "~{bam.file}" "~{bam.index}" \
            > ./sampleJsonUpdated.json  
    >>>

    output {
        SampleDescriptor sampleUpdated = read_json('./sampleJsonUpdated.json')
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
        Int timeMinutes = 5 + ceil(size(fileList, "G")) * 30
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
        timeMinutes: timeMinutes
    }
}
task ZipFiles {
    input {
        Array[File] fileList
        Array[String] optionalFileList = ["test.txt"]
        String outputPrefix
        Int memory = 512
        Boolean flattenArchive = true
        Int timeMinutes = 5 + ceil(size(flatten([fileList,select_all(optionalFileList)]), "G")) * 30
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
        timeMinutes: timeMinutes
    }
}

task CheckModules {
    input {
        Array[String] modules
        Int memory = 512
    }
    
    #WIP: does not localise optional files...

    command {
        set -e -o pipefail
        echo "# env"
        printenv > printenv.log
        echo "#modules"
        cat ~{write_lines(modules)} | 
        ( while read FILE; do 
            echo "Testing for the availability of $FILE"
            ml --quiet av $FILE;
            (ml $FILE || ml --ignore-cache load $FILE)
        done ) > modulelist.log

        
    }

    output {
        File moduleLog = "modulelist.log"
        File printenvLog = "printenv.log"
    }

    runtime {
        memory: memory
        timeMinutes: 20
    }
}