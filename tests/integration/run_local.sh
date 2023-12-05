(
    set -ex
    bash tests/run_project.sh \
        -i $PWD/tests/integration/json/fastqToVariants/inputs_local.json \
        -s $PWD/tests/integration/json/fastqToVariants/samples.json \
        -w $PWD \
        -r $PWD/tests/runs/integration_local \
        -f $PWD/tests/data/raw/fastq/ \
        -d /groups/umcg-pmb/tmp01//apps
)