# Bestie

WDL based hts-analysis for slurm cluster with enviroment modules. More scalable then ever before setting discrete cpu, memory, runtimes and disk space based on input files for maximum* sheduling efficiency.

> *: Always WIP due to lfs stability and edge cases.

### FastqToBam workflow

Shown in the example below

[![](https://mermaid.ink/img/pako:eNp9Ve1u4jgUfRXL0kisBqoWaIH82JUJhbYztIXQdlrTHy65hWgSO-M422Fp_84DzCPuk-y1E8rHrAYkRHzOPff4-l5nRWcqBOrR51i9zBZCGzLpTWXxJfhh_FymuSF9FkxGGflIhmBEKIx4JLXan6TLfQ3CAMmWSRzJr9ljEdZ1qM-ntC8yM_JJRYsXgtQw-2NKH7cz-I7as1Sbg_z74yfJu2xIKs_zp0hhyutoJnS4H9dzcac8mAljQBP1N_64wAz3EWY75FNH7vNAJBNlLX0rjfYdMOATHSUDESsNO8AZ2vJzu-PUkIpKTaSkiPetDBz3nAcQw8wQg1oJhEXNSrmzgrIdde6WLnj3jhEWR3OZgDS2wkJ_tXaSSM7L6AtH_cQDuzOCe9zd3ScHf-ZD0HNwcCnTy9MYi2dgfS6fHfOSBwpPGok7Mh8-kJvhOXnSQs4WxdKl41-t7DpI8RRD-NdbAV1Z6PUesldyzRGviReh4f_TXjuZkaWRkf-4LXCpXsmY-0qaSOaw78cqwfcdSyOnFaxKaM9VsHE1wbPbqkhFKhJCmKfF6VnyxEnd8IFWeTq2zdldosUSvXHoLfdFHK-NoM8MZJZnxNFL5q1j3mE-V9fKtxz0UooENqnuHOXLrw34xQH3v-kzy7p3rAdk2W7BOXrvl0ronG2owaase8UsipEKs8iKtQenyhi_1lDrjoJx2RIWG5fYngaO8hM845gQyy_vCFZwu2jPjnpaqu2PCYZvBRVXBPNX41y69fURMn9zhqzHuyKDMcxwx9gFRunSHyvmn51ylqbx0gqskWLYWX9tR2Vmy89WClskNuD9CGv9yzSw4gZgg53FYtTZGb_KDd6LmWfjcNowkbsbdTSzXUGrdK6jkHpG51ClCehE2Ee6sjpTahaQwJR6-DfEgZnSqXzDmFTIB6WSdRi25XxBvWcRZ_iUp3jvQi8Scy2S91UNMgTtq1wa6h21WsdOhXor-p169Xbr4LBdbzWbh-1O_bDZaFXpknq146PmQaPRaNWb9XbnuN5uvFXpPy5x4-Ck3ezUTxqd5lG9fdI5alYphBHWfVi8KNz74u0_7JraAg?type=png)](https://mermaid.ai/live/edit#pako:eNp9Ve1u4jgUfRXL0kisBqoWaIH82JUJhbYztIXQdlrTHy65hWgSO-M422Fp_84DzCPuk-y1E8rHrAYkRHzOPff4-l5nRWcqBOrR51i9zBZCGzLpTWXxJfhh_FymuSF9FkxGGflIhmBEKIx4JLXan6TLfQ3CAMmWSRzJr9ljEdZ1qM-ntC8yM_JJRYsXgtQw-2NKH7cz-I7as1Sbg_z74yfJu2xIKs_zp0hhyutoJnS4H9dzcac8mAljQBP1N_64wAz3EWY75FNH7vNAJBNlLX0rjfYdMOATHSUDESsNO8AZ2vJzu-PUkIpKTaSkiPetDBz3nAcQw8wQg1oJhEXNSrmzgrIdde6WLnj3jhEWR3OZgDS2wkJ_tXaSSM7L6AtH_cQDuzOCe9zd3ScHf-ZD0HNwcCnTy9MYi2dgfS6fHfOSBwpPGok7Mh8-kJvhOXnSQs4WxdKl41-t7DpI8RRD-NdbAV1Z6PUesldyzRGviReh4f_TXjuZkaWRkf-4LXCpXsmY-0qaSOaw78cqwfcdSyOnFaxKaM9VsHE1wbPbqkhFKhJCmKfF6VnyxEnd8IFWeTq2zdldosUSvXHoLfdFHK-NoM8MZJZnxNFL5q1j3mE-V9fKtxz0UooENqnuHOXLrw34xQH3v-kzy7p3rAdk2W7BOXrvl0ronG2owaase8UsipEKs8iKtQenyhi_1lDrjoJx2RIWG5fYngaO8hM845gQyy_vCFZwu2jPjnpaqu2PCYZvBRVXBPNX41y69fURMn9zhqzHuyKDMcxwx9gFRunSHyvmn51ylqbx0gqskWLYWX9tR2Vmy89WClskNuD9CGv9yzSw4gZgg53FYtTZGb_KDd6LmWfjcNowkbsbdTSzXUGrdK6jkHpG51ClCehE2Ee6sjpTahaQwJR6-DfEgZnSqXzDmFTIB6WSdRi25XxBvWcRZ_iUp3jvQi8Scy2S91UNMgTtq1wa6h21WsdOhXor-p169Xbr4LBdbzWbh-1O_bDZaFXpknq146PmQaPRaNWb9XbnuN5uvFXpPy5x4-Ck3ezUTxqd5lG9fdI5alYphBHWfVi8KNz74u0_7JraAg)


### BamToVariantsWorkFlow

Two workflows a generic germline workflow and a somatic workflow. The somatic workflow is subdivided into freebayes, MuTect2 (gatk4) and LoFreq.

Germline 

[![](https://mermaid.ink/img/pako:eNpVUluPojAU_itNk0lms2pAvCAPmyBWZINKkJlkl_rQlarsQEtK2d1Z43-fUoe5PJCcc_huPe0FHnhGoQOPBf97OBMhQbLADLO7OxBsoocEMzcNWNVIB-xIWRXU4-yYn8BXENMjFZQdqKoTIk5UBkxS8YcU9f5VIYpR5MYIzN31TgmBfv8bmF92ByIVENRar75iNtd_vDTM2VML3nejRTcKOsmVG4Xb5EeEPDcMUQwi9e3cdRQiFSPYJCh-dEPMPE1Hb15SBwR5l1CZLm6QVhfpcpliuCJVweVzRT1SFIp37z96S1CqHX3BsMuwRrGPPjhjttQCfurx8lfOaEuqQfV2xs9EFfyVqJbia-ZKWX_i3iv7bkEfnH20aQ8fbHzMVpoZpD5lOrEmdsDtQ6IvL9Cg7-kyZ6QA20aqm6wdjFkfRIL_pgcJWp4e3DCq28MePIk8g44UDe3BkoqStC28YAYAhvJMS4qho8qMiCcMMbsqTkXYT87LjiZ4czpD56jWrbqmyoiki5ycBHmHUJZR4fGGSeiY9mSkRaBzgf-gY9mzwcy0xqZhTe3h2LR68FmhLHNgTQ1jYpjDkWWMTfvag_-1rTGYju3J2JoNDWs0MUZDuwdplksu1rdHrt_69QU6eeoF?type=png)](https://mermaid.ai/live/edit#pako:eNpVUluPojAU_itNk0lms2pAvCAPmyBWZINKkJlkl_rQlarsQEtK2d1Z43-fUoe5PJCcc_huPe0FHnhGoQOPBf97OBMhQbLADLO7OxBsoocEMzcNWNVIB-xIWRXU4-yYn8BXENMjFZQdqKoTIk5UBkxS8YcU9f5VIYpR5MYIzN31TgmBfv8bmF92ByIVENRar75iNtd_vDTM2VML3nejRTcKOsmVG4Xb5EeEPDcMUQwi9e3cdRQiFSPYJCh-dEPMPE1Hb15SBwR5l1CZLm6QVhfpcpliuCJVweVzRT1SFIp37z96S1CqHX3BsMuwRrGPPjhjttQCfurx8lfOaEuqQfV2xs9EFfyVqJbia-ZKWX_i3iv7bkEfnH20aQ8fbHzMVpoZpD5lOrEmdsDtQ6IvL9Cg7-kyZ6QA20aqm6wdjFkfRIL_pgcJWp4e3DCq28MePIk8g44UDe3BkoqStC28YAYAhvJMS4qho8qMiCcMMbsqTkXYT87LjiZ4czpD56jWrbqmyoiki5ycBHmHUJZR4fGGSeiY9mSkRaBzgf-gY9mzwcy0xqZhTe3h2LR68FmhLHNgTQ1jYpjDkWWMTfvag_-1rTGYju3J2JoNDWs0MUZDuwdplksu1rdHrt_69QU6eeoF)



### Goals:

- [x] is to get generic alignment working.
- [x] basic variant calling (haplotypecallerGvcf).
- [x] somatic variant calling(MuTect2,...).
- [ ] ichorCNA integration.
- [x] Variant annotation (vep)
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

| Name         | project website                                                                            | Article          |
| ------------ | ------------------------------------------------------------------------------------------ | ---------------- |
| GNU Parallel | [gnu.org](https://www.gnu.org/software/parallel/parallel_design.html) | [doi](https://doi.org/10.5281/zenodo.12789352) |
| Fastqc       | [bioinformatics.babraham.ac.uk](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | |
| BWA          | [github](https://github.com/lh3/bwa)                                                       | [preprint](http://arxiv.org/abs/1303.3997) |
| Picard-tools | [sourceforge](http://picard.sourceforge.net/)                                              | [instructions below faq](http://picard.sourceforge.net/) |
| GATK4 + MuTect2 toolkit | [project home](http://www.broadinstitute.org/gatk/)                                        | [instructions here](https://www.broadinstitute.org/gatk/about/citing-gatk) |
| SAMtools     | [project home](http://www.htslib.org/)                                                     | [pubmed](http://www.ncbi.nlm.nih.gov/pubmed/19505943)|
| HTSeq		| [project home](http://www-huber.embl.de/users/anders/HTSeq/doc/index.html)		    | [pubmed](http://www.ncbi.nlm.nih.gov/pubmed/25260700) |
| Cutadapt     | [github](https://github.com/marcelm/cutadapt/) | [doi](http://dx.doi.org/10.14806/ej.17.1.200) | 
| bcftools     | [github](https://samtools.github.io/bcftools/) | [pubmed](https://pubmed.ncbi.nlm.nih.gov/33590861) | 
| fgbio        | [github](https://fulcrumgenomics.github.io/fgbio/tools/latest/) |  | 
| freebayes    | [github](https://github.com/freebayes/freebayes) | [arxiv.org](http://arxiv.org/abs/1207.3907) |
| IchorCNA     | [github](https://github.com/broadinstitute/ichorCNA/) | [doi](https://doi.org/10.1038/s41467-017-00965-y) | 
| LoFreq       | [github](https://github.com/csb5/lofreq/) | [pubmed](http://www.ncbi.nlm.nih.gov/pubmed/23066108) | 
| MarkTrimming |[github](https://github.com/mmterpstra/marktrimming) |  |
| MultiQC      | [seqera](https://docs.seqera.io/multiqc) | [doi](http://dx.doi.org/10.1093/bioinformatics/btw354) | 
| VEP          |[github](https://github.com/Ensembl/ensembl-vep) | [doi](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0974-4) | 
| TrimGalore   |[github](https://github.com/FelixKrueger/TrimGalore) |  | 
| GATK Bundle   | [ human reference ](http://gatkforums.broadinstitute.org/discussion/1213/what-s-in-the-resource-bundle-and-how-can-i-get-it) |
| Ensembl       | [reference/gtf dowload](http://www.ensembl.org/info/data/ftp/index.html) |
| UCSC Tools    | [ format conversion/additional tools ](http://hgdownload.soe.ucsc.edu/admin/exe/) |
