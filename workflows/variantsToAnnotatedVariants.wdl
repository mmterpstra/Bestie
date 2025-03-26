version 1.0

import "../structs.wdl" as structs
import "../tasks/common.wdl" as common
import "../tasks/gatk.wdl" as gatk

workflow VariantsToAnnotatedVariants {
    input {
        String picardModule = "picard/2.26.10-Java-8-LTS"
        String gatkModule = "GATK/4.2.4.1-Java-8-LTS"
        #String samtoolsModule = "SAMtools/1.15.1-GCC-11.3.0"
        Reference reference
        File funcotatoTar
        IndexedFile variants
        Boolean somatic = true
        File? sampleJson
        SampleConfig? sampleConfigIn 
        File targetIntervalList
        Int targetScatter
    }
    variantsFile = variants.file
    call gatk.FuncotateVcfs as funcotateVariants {
        input:
            inputVcf= variants,
            inputVcfsFile=variantsFile,
            
    }
    
}