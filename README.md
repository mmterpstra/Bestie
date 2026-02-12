# Bestie

WDL based hts-analysis for slurm cluster with enviroment modules. More scalable then ever before setting discrete cpu, memory, runtimes and disk space based on input files for maximum* sheduling efficiency.

> *: Always WIP due to lfs stability and edge cases.

### Goals:

- [x] is to get generic alignment working.
- [x] basic variant calling (haplotypecallerGvcf).
- [x] somatic variant calling(MuTect2,...).
- [ ] ichorCNA integration.
- [ ] Variant annotation
- [ ] Functional filtering of vcfs
- [x] Stability over > 10 samples
- [ ] Rework dir structure to be more in line with warp/other public resources
- [ ] End to end easybuild install

Relevant example resources for writing wdl files :

[WARP](https://github.com/broadinstitute/warp)
[BioWDL](https://github.com/biowdl/)

### How to use

Run the cromwell validation tool womtool to validate the input and generate a template input file.
for sampleJson see ./tests/data/raw/fastq/samples.json. See below for the most simple run example 

```
 java -Xmx8g -Dconfig.file=./path/to/cromwell.conf -jar ./path/to/cromwell.jar run Bestie.wdl -i inputs_integration.json
```

#### prepping a samplesheet/json for input

This first part makes a sample tsv with one readgroup per line that can be edited as needed 

```
ls /path/to/raw/fastq/*_R1.fastq.gz | perl -wne 'BEGIN{print "fq1,fq2,sampleName\n"};chomp;print $_;s/_R1./_R2./g; print ",$_"; s/.*/([\w\d]*)_.*/$1/g; print ",$_" ;print "\n"'| perl SampleSheetTool.pl reformatmin /dev/stdin > samplesheet.csv
```
Converts the samplesheet to json file merging readgroups to samples as needed (some text alignment issues)

```
perl SamplesheetTool jsondump samplesheet.csv > samplesheet.json
```

The next command sets up the environment in the /path/to/workflow/output_folder/ copying all the needed files and linking all the needed folders raw and runs `cromwell.jar` to execute the workflow.

- This needs '-d /path/to/folder/containing/data/' usually the 'apps/' folder for me.
- This results in outputs below  '/path/to/workflow/output_folder/' containing all the (fixupped) files needed for running the analysis.

```
fastq folders (example based on `tests/integration/run_local.sh`):
(
    set -ex
    bash tests/run_project.sh \
        -i $PWD/tests/integration/json/fastqToVariants/inputs_local.json \
        -s samplesheet.json \
        -w $PWD \
        -r /path/to/workflow/output_folder/ \
        -f /path/to/raw/fastq/ \
        -d /path/to/folder/containing/data/
)
```



### How to install

easybuild the required modules or use future wrapper module

### Used tools and databases 

| Tool/data | Site | Citation | Doi  |
| --------- | ---- | -------- | ---- | 
| GNU Parallel | https://www.gnu.org/software/parallel/parallel_design.html#citation-notice | https://www.gnu.org/software/parallel/parallel_design.html#citation-notice |  https://doi.org/10.5281/zenodo.12789352 |
| Name         | project website                                                                            | Article          |
| ------------ | ------------------------------------------------------------------------------------------ | ---------------- |
| Fastqc       | [bioinformatics.babraham.ac.uk](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | |
| BWA          | [github](https://github.com/lh3/bwa)                                                       | [preprint](http://arxiv.org/abs/1303.3997) |
| Picard-tools | [sourceforge](http://picard.sourceforge.net/)                                              | [instructions below faq](http://picard.sourceforge.net/) |
| GATK4 + MuTect2 toolkit | [project home](http://www.broadinstitute.org/gatk/)                                        | [instructions here](https://www.broadinstitute.org/gatk/about/citing-gatk) |
| SAMtools     | [project home](http://www.htslib.org/)                                                     | [pubmed](http://www.ncbi.nlm.nih.gov/pubmed/19505943)|
| HTSeq		| [project home](http://www-huber.embl.de/users/anders/HTSeq/doc/index.html)		    | [pubmed](http://www.ncbi.nlm.nih.gov/pubmed/25260700) |
| Cutadapt     |  |  | 
| bcftools     |  |  | 
| fgbio        |  |  | 
| freebayes    |  |  |
| IchorCNA     |  |  | 
| LoFreq       |  |  | 
| MarkTrimming |  |  |
| MultiQC      |  |  | 
| VEP          |  |  | 
| TrimGalore   |  |  | 
| GATK Bundle   | [ human reference ](http://gatkforums.broadinstitute.org/discussion/1213/what-s-in-the-resource-bundle-and-how-can-i-get-it) |
| Ensembl       | [reference/gtf dowload](http://www.ensembl.org/info/data/ftp/index.html) |
| UCSC Tools    | [ format conversion/additional tools ](http://hgdownload.soe.ucsc.edu/admin/exe/) |
