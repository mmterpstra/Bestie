#execute this from the root of the project `bash src/run.sh`
#(cd wdl-bamstats/&& ml cromwell/56-Java-11 && java -jar $EBROOTCROMWELL/cromwell.jar run  Dockstore.wdl -i test.wdl.json --workflow-root ./ )
#when running with slurm as a backend
#(cd wdl-bamstats/ && ml cromwell/56-Java-11 && java -Dconfig.file=../siteconfig/peregrine/slurm.conf -jar $EBROOTCROMWELL/cromwell.jar run Dockstore.wdl -i test.wdl.json --workflow-root ./ )


(
    #init env
    WORKFLOWROOT="$(pwd)"
    #echo $WORKFLOWROOT
    RUNROOT='./tests/runs/'
    mkdir -p $RUNROOT/raw
    ln -s $(pwd)/tests/data/raw/fastq $RUNROOT/raw
    ln -s $(pwd)/tests/data/ref $RUNROOT/ref
    if [ -e $RUNROOT/inputs.json ] ; then
        rm -v $RUNROOT/inputs.json
    fi
    ln -s $(pwd)/inputs_integration.json $RUNROOT/inputs.json
    #start workflow
    set -x
    ml cromwell/56-Java-11 && \
    java -Xmx8g -Dconfig.file=$WORKFLOWROOT/site/gearshift/cromwell.conf -jar $EBROOTCROMWELL/womtool.jar validate $WORKFLOWROOT/Bestie.wdl -i $RUNROOT/inputs.json 
    cd $RUNROOT
    java -Xmx8g -Dconfig.file=$WORKFLOWROOT/site/gearshift/cromwell.conf -jar $EBROOTCROMWELL/cromwell.jar run $WORKFLOWROOT/Bestie.wdl -i inputs.json --workflow-root $WORKFLOWROOT
)
