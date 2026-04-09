(
    set -ex

    bash tests/run_project.sh \
        -i $PWD/tests/integration/json/fastqToVariants/inputs_local_twist_umi_consensus.json \
        -s $PWD/tests/integration/json/fastqToVariants/samples_umi.json \
        -w $PWD \
        -r $PWD/tests/runs/integration_local \
        -f $PWD/tests/data/raw/fastq/ \
        -d $PWD/tests/
)
