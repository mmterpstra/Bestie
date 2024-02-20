#execute this from the root of the project `bash src/run.sh`
#(cd wdl-bamstats/&& ml cromwell/56-Java-11 && java -jar $EBROOTCROMWELL/cromwell.jar run  Dockstore.wdl -i test.wdl.json --workflow-root ./ )
#when running with slurm as a backend
#(cd wdl-bamstats/ && ml cromwell/56-Java-11 && java -Dconfig.file=../siteconfig/peregrine/slurm.conf -jar $EBROOTCROMWELL/cromwell.jar run Dockstore.wdl -i test.wdl.json --workflow-root ./ )

#smart defaults
PIPELINE="Bestie.wdl"
WORKFLOWROOT="$PWD"
CONFIG=$(ls $WORKFLOWROOT/site/*/cromwell.conf) 

while getopts i:s:w:r:f:d:p:c: flag
do
    case "${flag}" in
        i) INPUTJSON=${OPTARG};;
        s) SAMPLEJSON=${OPTARG};;
        w) WORKFLOWROOT=${OPTARG};;
        r) RUNROOT=${OPTARG};;
        f) FASTQRAWDIRS=${OPTARG};;
        d) DATADIR=${OPTARG};;
        p) PIPELINE=${OPTARG};;
        c) CONFIG=${OPTARG};;
    esac
done

ml purge
(
    #init env
    
    set -e
    #create limited acces raw folder
    mkdir -p $RUNROOT/raw

    #init data
    if [ ! -e "$RUNROOT/data" ]; then
        #ln -s /groups/umcg-pmb/tmp01//apps/data $RUNROOT/
        ln -s $DATADIR/data $RUNROOT/
        
    fi

    #archive inputs.json
    if [ ! -e "$RUNROOT/inputs.json" ]; then
        cp $INPUTJSON $RUNROOT/inputs.json
        cp $INPUTJSON $SAMPLESHEETFOLDER/${HOSTNAME}_$(basename $INPUTJSON)
    elif [ -e "$RUNROOT/inputs.json" ]; then
        mv $RUNROOT/inputs.json{,$(date --iso-8601=min | perl -wpe 's/[\+\:]/_/g').bak}
        cp $INPUTJSON $RUNROOT/inputs.json
        echo "## "$(date --iso-8601=min | perl -wpe 's/[\+\:]/_/g')"## $PWD $0 $@" >> $SAMPLESHEETFOLDER/${HOSTNAME}_bestie.log
        
        #echo $SAMPLESHEETFOLDER/${HOSTNAME}_$(basename $INPUTJSON){,$(date --iso-8601=min | perl -wpe 's/[\+\:]/_/g').bak}
        mv $SAMPLESHEETFOLDER/${HOSTNAME}_$(basename $INPUTJSON){,$(date --iso-8601=min | perl -wpe 's/[\+\:]/_/g').bak}
        cp $INPUTJSON $SAMPLESHEETFOLDER/${HOSTNAME}_$(basename $INPUTJSON)
    fi

    #archive sample.json
    if [ ! -e "$RUNROOT/sample.json" ]; then
        cp $SAMPLEJSON $RUNROOT/sample.json
        cp $INPUTJSON $SAMPLESHEETFOLDER/${HOSTNAME}_$(basename $SAMPLEJSON)
    elif [ -e "$RUNROOT/sample.json" ]; then
        mv $RUNROOT/sample.json{,$(date --iso-8601=min | perl -wpe 's/[\+\:]/_/g').bak}
        cp $SAMPLEJSON $RUNROOT/sample.json
        mv $SAMPLESHEETFOLDER/${HOSTNAME}_$(basename $SAMPLEJSON){,$(date --iso-8601=min | perl -wpe 's/[\+\:]/_/g').bak}
        cp $INPUTJSON $SAMPLESHEETFOLDER/${HOSTNAME}_$(basename $SAMPLEJSON)
    fi
    
    for FASTQRAWDIR in $(echo $FASTQRAWDIRS | tr , \ ); do
        if [ ! -e "$RUNROOT/raw/$(basename $FASTQRAWDIR)" ]; then
            ln -s $FASTQRAWDIR $RUNROOT/raw
        else 
            unlink $RUNROOT/raw/$(basename $FASTQRAWDIR)
            ln -s $FASTQRAWDIR $RUNROOT/raw
        fi
        perl -i$(basename $FASTQRAWDIR).bak -wpe 's?'$FASTQRAWDIR'?'$RUNROOT/raw/$(basename $FASTQRAWDIR)/'?g' $RUNROOT/sample.json
    done
    #fix samplejson link?
    #ls -alh $RUNROOT
    perl -i.bak -wpe 's?.sampleJson": ".*",?.sampleJson": "'$RUNROOT/sample.json'",?g' $RUNROOT/inputs.json

    #start workflow
    ml cromwell/79-Java-11|| ml cromwell
    
    set -x
    java -Xmx8g -Dconfig.file=$CONFIG \
        -jar $EBROOTCROMWELL/womtool.jar validate \
        $WORKFLOWROOT/$PIPELINE \
        -i $RUNROOT/inputs.json 
    cd $RUNROOT
    
    echo "bash ${SCRIPTCALL}," $(git log -1 || echo "$(pwd);$(date)" ) | head -n 1  >> "$RUNROOT/nohup_"$(date --rfc-3339=date)".out"

    nohup java -Xmx8g -Dconfig.file=$CONFIG \
        -jar $EBROOTCROMWELL/cromwell.jar run \
        $WORKFLOWROOT/$PIPELINE \
        -i $RUNROOT/inputs.json \
        --workflow-root $WORKFLOWROOT &>> "$RUNROOT/nohup_"$(date --rfc-3339=date)".out"
    echo "watch log in $RUNROOT/nohup_"$(date --rfc-3339=date)".out"
)