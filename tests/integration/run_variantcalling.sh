(
    #this should test "workflows/bamsToVariants.wdl" without runnning the fastqToBam workflow first
    bash tests/run_project.sh \
        -i $PWD/tests/integration/json/bamToVariants/inputs.json \
        -s $PWD/tests/integration/json/bamToVariants/sampleswithbam.json \
        -w $PWD \
        -r $PWD/tests/runs/variantstoBam \
        -f $PWD/tests/data/raw/fastq/ \
        -d $PWD/tests/ \
        -p workflows/bamsToVariants.wdl
)