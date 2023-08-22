# Bestie

wdl based hts-analysis for slurm cluster with enviroment modules. 

### Goals:

- [x] is to get generic alignment working.
- [x] basic variant calling (haplotypecallerGvcf).
- [ ] somatic variant calling(MuTect2,...).
- [ ] ichorCNA integration.
- [ ] Variant annotation
- [ ] Functional filtering of vcfs
- [ ] Stability over > 10 samples
- [ ] Rework dir structure to be more in line with warp/other public resources

Relevant example resources for writing wdl files :

[WARP](https://github.com/broadinstitute/warp)
[BioWDL](https://github.com/biowdl/)

### How to use

Run the cromwell validation tool womtool to validate the input and generate a template input file.
for sampleJson see ./tests/data/raw/fastq/samples.json. See below for the most simple run example 

```
 java -Xmx8g -Dconfig.file=./path/to/cromwell.conf -jar ./path/to/cromwell.jar run Bestie.wdl -i inputs_integration.json
```

### How to install

easybuild the required modules or use future wrapper module
