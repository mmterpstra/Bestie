(
    #this should test "workflows/bamsToVariants.wdl" without runnning the fastqToBam workflow first
    bash tests/run_project.sh \
        -i $PWD/tests/integration/json/bamToVariants/inputs_local.json \
        -s $PWD/tests/integration/json/bamToVariants/sampleswithbam_local.json \
        -w $PWD \
        -r $PWD/tests/runs/bamToVariants_local \
        -f $PWD/tests/data/raw/fastq/,$PWD/tests/data/raw/bam/ \
        -d /groups/umcg-pmb/tmp01//apps \
        -p workflows/bamsToVariants.wdl
)