(
    bash tests/run_project.sh \
        -i $PWD/tests/integration/json/fastqToVariants/inputs.json \
        -s $PWD/tests/integration/json/fastqToVariants/samples.json \
        -w $PWD \
        -r $PWD/tests/runs/integration \
        -f $PWD/tests/data/raw/fastq/ \
        -d $PWD/tests/
)