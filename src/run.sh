#execute this from the root of the project `bash src/run.sh`
(cd wdl-bamstats/ &&  java -jar ../cromwell-84.jar run  Dockstore.wdl -i test.wdl.json --workflow-root ./ )
#when running with slurm as a backend
#(cd wdl-bamstats/ &&  java -Dconfig.file=../peregrine/slurm.conf -jar ../cromwell-84.jar  Dockstore.wdl -i test.wdl.json --workflow-root ./ )