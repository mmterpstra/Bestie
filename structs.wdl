version 1.0

struct IndexedFile {
    File file
    File index
}
struct ReadGroup {
    #should be unique for multiple samples??
    String identifier
    File fastq1
    File? fastq2
    File? fastqUmi
    
    #docs for read structures: https://fulcrumgenomics.github.io/fgbio/tools/1.5.1/ExtractUmisFromBam.html
    # Read structures are made up of <number><operator> pairs much like the CIGAR string in BAM files. Four kinds of operators are recognized:
    #    T identifies a template read
    #    B identifies a sample barcode read
    #    M identifies a unique molecular index read
    #    S identifies a set of bases that should be skipped or ignored

    String? readStructure1
    #default +T
    String? readStructure2
    #default +M 
    String? readStructureUmi 

    #usually assume library is approximated by sample_barcode1(+barcode2) if not here is your parameter to overwrite
    String? library
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
    Boolean? runTwistUmi
    IndexedFile? alignedReads
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
    File? altFile
}

struct Reference {
    File fasta
    File dict
    File fai
}