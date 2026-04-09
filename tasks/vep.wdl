version 1.0

import "../structs.wdl" as structs

task AnnotateTmpStorage {
    input {
        IndexedFile inputIdxVcf
        String outputBasename
        Reference reference
        File vepCacheTar
        #File vepCacheTar

        String vepModule = "VEP"
        Int memoryGb = 8
        Int threads = 2
        Int timeMinutes = 25 + ceil(size(inputIdxVcf.file, "G")) * 120
        Int disk = ceil(size(vepCacheTar, "M")*1.1)
        Boolean byReadGroup = false
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
        set -eo pipefail
        echo "$(date) ## INFO ## Unpacking cache tar"
        (cd $TMPDIR && tar -xf ~{vepCacheTar}) 
        echo "$(date) ## INFO ## module load"
        ml ~{vepModule}
        
        ls $TMPDIR/113
        ls  -alh $TMPDIR/*
        echo "$(date) ## INFO ## Running vep"
        vep \
        --input_file ~{inputIdxVcf.file} \
        --species homo_sapiens --refseq \
        --cache --offline \
        --output_file ~{outputBasename}.vcf \
        --dir_cache $TMPDIR/113 \
        --use_given_ref \
        --fasta ~{reference.fasta} \
        --force_overwrite \
        --everything \
        --hgvs \
        --dir_plugins $TMPDIR/113/Plugins \
        --exclude_predicted \
        --symbol \
        --flag_pick_allele \
        --total_length \
        --shift_3prime "1" \
        --allele_number \
        --numbers \
        --dont_skip \
        --allow_non_variant \
        --buffer_size 1000 \
        --fork 2 \
        --plugin "SpliceAI,snv=$TMPDIR/113/Plugins/data/hg38/spliceai_scores.raw.snv.ensembl_mane.grch38.107.vcf.gz,indel=$TMPDIR/113/Plugins/data/hg38/spliceai_scores.raw.snv.ensembl_mane.grch38.107.vcf.gz" \
        --plugin AlphaMissense,file=$TMPDIR/113/Plugins/data/hg38/AlphaMissense_hg38.tsv.gz \
        --plugin CADD,snv=$TMPDIR/113/Plugins/data/hg38/whole_genome_SNVs.tsv.gz,indels=$TMPDIR/113/Plugins/data/hg38/gnomad.genomes.r4.0.indel.tsv.gz \
        --vcf
        #        --plugin "PolyPhen_SIFT,dir=$TMPDIR/113/Plugins/data/hg38/"


        echo "$(date) ## INFO ## creating bgzipped and tabixed output"
        bgzip ~{outputBasename}.vcf
        tabix -p vcf ~{outputBasename}.vcf.gz
        echo "$(date) ## INFO ## Cleanup of tmpdir"
        rm -r $TMPDIR/*
        echo "$(date) ## INFO ## done"
    }
    
    output {
        File vcf = outputBasename + ".vcf.gz"
        File idx = outputBasename + ".vcf.gz.tbi"
        IndexedFile idxVcf = { 
          "file" : outputBasename + ".vcf.gz",
          "index" : outputBasename + ".vcf.gz.tbi"
        }
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
        cpus: threads
    }
}

task InstallVepRefseqCache {
    input {
        File inputBam
        File targetIntervalList
        String outputMetricsBasename
        Reference reference
        String vepModule = "VEP"
        Int? memoryGb = "4"
        Int timeMinutes = 1 + ceil(size(inputBam, "G")) * 120
        Int disk = ceil(size(inputBam, "M")*2.1)
        Boolean byReadGroup = false
    }
    #https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl#L96
    command {
        #This aint ment to be ran just notes on how to install 
        ml ~{vepModule}
        mkdir -p ./apps/data/vep
        export VEP_DIR="./apps/data/vep/113/"
        export PLUGINS="$VEP_DIR/plugins"
        export PLUGINDATA="$Plugins/data/hg38"
        mkdir -p $PLUGINDATA
        export CACHEDIR="$VEP_DIR/cache"
        mkdir -p $CACHEDIR
        (   
            ml VEP 
            INSTALL.pl -a p cf -s homo_sapiens_refseq -y GRCh38 \
            --PLUGINS AlphaMissense,SpliceAI,Grantham,gnomAD,PolyPhen_SIFT,gnomADc,SpliceRegion,SpliceVault,CADD,CAPICE \
            --CACHE_DIR $CACHEDIR
        wget https://github.com/molgenis/vip/raw/refs/heads/main/resources/vep/plugins/Grantham.pm -O $PLUGINS/Plugins/Grantham.pm
        wget https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz \
          -O $PLUGINDATA/whole_genome_SNVs.tsv.gz 
        wget https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz.tbi \
          -O $PLUGINDATA/whole_genome_SNVs.tsv.gz.tbi
        wget https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz \
            -O $PLUGINDATA/gnomad.genomes.r4.0.indel.tsv.gz
        wget https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz \
            -O $PLUGINDATA/gnomad.genomes.r4.0.indel.tsv.gz.tbi
        wget https://storage.cloud.google.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz -O $PLUGINDATA/AlphaMissense_hg38.tsv.gz
        tabix -s 1 -b 2 -e 2 -f -S 1 $PLUGINDATA/AlphaMissense_hg38.tsv.gz
        wget https://storage.cloud.google.com/dm_alphamissense/AlphaMissense_gene_hg38.tsv.gz -O $PLUGINDATA/AlphaMissense_gene_hg38.tsv.gz
        #SpliceAI plugin needs the following files: source: https://github.com/broadinstitute/SpliceAI-lookup/blob/master/README.md#local-install  
        wget https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_plugins/spliceai_scores.raw.snv.ensembl_mane.grch38.107.vcf.gz \
            -O  $PLUGINDATA/spliceai_scores.raw.snv.ensembl_mane.grch38.107.vcf.gz
        wget https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_plugins/spliceai_scores.raw.snv.ensembl_mane.grch38.107.vcf.gz.tbi \
            -O  $PLUGINDATA/spliceai_scores.raw.snv.ensembl_mane.grch38.107.vcf.gz.tbi

        #should work like this:
        #(ml VEP && vep \
        #--input_file in.vcf.gz \
        #--species homo_sapiens --refseq \
        #--cache --offline \
        #--output_file out.vcf.gz \
        #--dir_cache .../apps/data//VEP/113/ \
        #--fasta .../genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
        #--force_overwrite \
        #-use_given_ref \
        #--everything 
        # --dir_plugins $PLUGINS ) 
        (ml VEP && vep --input_file scat_8.vcf.gz \
        --species homo_sapiens --refseq --cache --offline \
         --output_file test_8.vcf.gz --dir_cache ./apps/data//vep/113/  --fasta ./Homo_sapiens_assembly38.fasta --force_overwrite --use_given_ref \
         --everything --dir_plugins ./apps/data//vep/113/plugins/Plugins --exclude_predicted \
         --symbol --flag_pick_allele --total_length --shift_3prime "1" --allele_number --numbers \
         --dont_skip --allow_non_variant --buffer_size 1000 --fork 2 \
         --plugin PolyPhen_SIFT,dir=/some/specific/directory
         --plugin "SpliceAI,snv=./apps/data//vep/113/Plugins/data/hg38/spliceai_scores.raw.snv.ensembl_mane.grch38.107.vcf.gz,indel=./apps/data//vep/113/Plugins/data/hg38/spliceai_scores.raw.snv.ensembl_mane.grch38.107.vcf.gz,cutoff=0.5" \
         --plugin AlphaMissense,file=./apps/data//vep/113/Plugins/data/hg38/AlphaMissense_hg38.tsv.gz \
         --plugin CADD,snv=./apps/data/vep/113/Plugins/data/hg38/whole_genome_SNVs.tsv.gz,indels=./apps/data//vep/113/Plugins/data/hg38/gnomad.genomes.r4.0.indel.tsv.gz \
         --vcf  )
    }
    
    output {
        File hsMetrics = outputMetricsBasename + ".hs_metrics"
    }

    runtime {
        memory: select_first([memoryGb * 1024,4*1024])
        timeMinutes: timeMinutes
        disk: disk
    }
}