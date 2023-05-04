version 1.0

struct ReadGroup {
    #should be unique for multiple samples??
    String identifier
    File fastq1
    File? fastq2
    File? fastq_umi
    #[ATCGN]{6+}
    String barcode1
    String? barcode2
    #001 or something
    String? run
    #HXXXXXXXX or something 
    String? flowcell
    #L000
    String? lane
    #YYYYMMDD
    String? date
    #hiseq/miseq etc
    String? sequencer
    #illumina
    String? platform
    #
    String? platform_model
    #
    String? sequencing_center
}

struct SampleDescriptor {
    String name
    String? threeLetterName
    String? control
    String? gender
    Array[ReadGroup] readgroups
}


struct SampleConfig {
    Array[SampleDescriptor] samples
}

struct BwaIndex {
    File fastaFile
    File ambFile
    File annFile
    File bwtFile
    File pacFile
    File saFile
    #placeholder for optional alt stuff
    File? altFiles
}

struct Reference {
    File fasta
    File dict
    File fai
}