# Bestie

WDL based hts-analysis for slurm cluster with enviroment modules. More scalable then ever before setting discrete cpu, memory, runtimes and disk space based on input files for maximum* sheduling efficiency.

> *: Always WIP due to lfs served over nfs stability and edge cases.

### FastqToBam workflow

Shown in the example below. In short an pipeline for analysising reads from illumina sequecers using twist udi adapters in a cfDNA setting. Using parallelisation if possible. One of the features is the adapter trimming and alignment step substep is parallelised in batches smaller then 50m reads. 

<!---
    The workflow should look like:
        flowchart TD

        A[Input FASTQs + Metadata] --\> B[Create symlinks]
        B --\> C["FastQC (raw reads)"]

        B --\> D["FASTQ → uBAM (fgbio + Picard)"]

        D --\> E[Scatter over uBAM shards]
        E --\> F[SamToFastq]
        F --\> G[TrimGalore]
        F --\> H["Cutadapt (optional)"]


        G --\> I[Select trimmed FASTQ]
        H --\> I

        F --\> J[BWA Alignment + MarkTrimming]
        I --\> J
        J --\> K[Shard BAMs]

        K --\> L[Merge BAMs + MarkDuplicates]
        L --\> N[Sort BAM]

        %% UMI branch
        K --\> O{UMI enabled?}
        O --\>|Yes| P[UMI-aware MarkDuplicates]
        P --\> Q[UMI QC]

        %% Duplex branch
        O --\>|Yes| S{Duplex enabled?}
        S --\>|Yes| T["Merge BAMs (no dedup)"]
        T --\> U[GroupReadsByUMI]
        U --\> V[Call Duplex Consensus Reads]
        V --\> W["Sort (queryname)"]
        W --\> X[SamToFastq]
        X --\> Y["Cutadapt (optional)"]
        X --\> Z["BWA realignment (duplex)"]
        Y --\> Z


        %% Merge paths
        Z --\> AA[Pre-BQSR BAM]
        N --\> AA 
        P --\> AA

        %% QC before BQSR
        AA --\> AB["QC (pre-BQSR)"]


        %% BQSR
        AA --\> AC{Run BQSR?}
        AC --\>|Yes| AD[BaseRecalibrator]
        AD --\> AE[ApplyBQSR]
        AE --\> AF["QC (post-BQSR)"]
        AC --\>|No| AG[Final BAM]
        AE --\> AG

        AG --\> AH[Outputs: BAM + QC + Metrics]
        AB --\> AH
        AF --\> AH
        Q --\> AH
        C --\> AH

-->

[![](https://mermaid.ink/img/pako:eNp9Ve1O2zAUfRXLElKnFdQWaKE_NjkpLQwKlJQxcPlhGtNGJHbmOGNdx989wB5xT7JrO6FNmdZKVW_Oufce3w9niacy5LiLH2P5PJ0zpdG4NxHui-BD6IlIc436JBiPMvQeDblmIdPsHm1vf0Ae9RVnmqNskcSReMrunZtnUZ9OcJ9leuSjmmLPCKhh9m6C79czOGrPUE0O9OfXb5R7ZIhqj7OHSELKy2jKVLjp17N-RzSYMq25QvIb_FjHDM4RZhXykSX3acCSsTSSvhZC-xYY0LGKkgGLpeIV4Bhk-bk5capRTaY6koLFm1IGlntCAx7zqUYaYiU8dDUrwh07yrqXy_CJejcEkTiaiYQLbSrM1JORk0RiVnifOKozPlnjlAbmmAgOXD3qqYXP6JCrGbdwEbOXpzFUUvOySWeWeU4DCW0HYiXM1ha6Hp6gB8XEdL4e-WJpnnPBHmIefnxx0IWBft7y7Ce6pIBvs2em-L_TXtowI0NDI38zqaHz75W8a8GDZYFv5A9WlDG0bO3sNSFRyMM8dU0z5LEVcE0HSubplZlJbwFiCvTaop-pz-K4VONLkXGR5Rmy9IL52TJvIJ-tYO1rztVCsISvUt1Yype3c_fFArf_Ga8V6w5YZkhgfV7HpBZaZSvqraNuVNMVImV6nrlnd5ZGCL1UfNsbBVdF4w12XmBovU-EVALCKj_wR1gTZJyLO4I4pgc6zaqnRehS3Iamt47-8ioX9nnZT-KvGkp61GMZv-JTOD7MhZaqEEzcHUCOKEnTeGEClIhbeNIvJclMVzS9pjiXkGFA-xEU_s0SlGEGlYdu3ckxvcg13I1Z1_jBkkEiez-qaFqOCPEKbmH2K-aoYvmvlvniOp6pKMRdrXJexwlXCTMmXtqSYj3nCZ_gLvwNYcsmeCJewCdl4k7KpHSDCZ_NcfeRxRlYeQo3N-9FbKbYisJFyJUvc6Fxd393r22D4O4Sf8fd5l5np723d9DYP2wctBq7zVYdL4B2uNPpNNq7nf1mq9XoNHZf6viHTdvcOei02odtIHfahwedZh3zMIKWDd17xr5uXv4CRrzlLg?type=png)](https://mermaid.ai/live/edit#pako:eNp9Ve1O2zAUfRXLElKnFdQWaKE_NjkpLQwKlJQxcPlhGtNGJHbmOGNdx989wB5xT7JrO6FNmdZKVW_Oufce3w9niacy5LiLH2P5PJ0zpdG4NxHui-BD6IlIc436JBiPMvQeDblmIdPsHm1vf0Ae9RVnmqNskcSReMrunZtnUZ9OcJ9leuSjmmLPCKhh9m6C79czOGrPUE0O9OfXb5R7ZIhqj7OHSELKy2jKVLjp17N-RzSYMq25QvIb_FjHDM4RZhXykSX3acCSsTSSvhZC-xYY0LGKkgGLpeIV4Bhk-bk5capRTaY6koLFm1IGlntCAx7zqUYaYiU8dDUrwh07yrqXy_CJejcEkTiaiYQLbSrM1JORk0RiVnifOKozPlnjlAbmmAgOXD3qqYXP6JCrGbdwEbOXpzFUUvOySWeWeU4DCW0HYiXM1ha6Hp6gB8XEdL4e-WJpnnPBHmIefnxx0IWBft7y7Ce6pIBvs2em-L_TXtowI0NDI38zqaHz75W8a8GDZYFv5A9WlDG0bO3sNSFRyMM8dU0z5LEVcE0HSubplZlJbwFiCvTaop-pz-K4VONLkXGR5Rmy9IL52TJvIJ-tYO1rztVCsISvUt1Yype3c_fFArf_Ga8V6w5YZkhgfV7HpBZaZSvqraNuVNMVImV6nrlnd5ZGCL1UfNsbBVdF4w12XmBovU-EVALCKj_wR1gTZJyLO4I4pgc6zaqnRehS3Iamt47-8ioX9nnZT-KvGkp61GMZv-JTOD7MhZaqEEzcHUCOKEnTeGEClIhbeNIvJclMVzS9pjiXkGFA-xEU_s0SlGEGlYdu3ckxvcg13I1Z1_jBkkEiez-qaFqOCPEKbmH2K-aoYvmvlvniOp6pKMRdrXJexwlXCTMmXtqSYj3nCZ_gLvwNYcsmeCJewCdl4k7KpHSDCZ_NcfeRxRlYeQo3N-9FbKbYisJFyJUvc6Fxd393r22D4O4Sf8fd5l5np723d9DYP2wctBq7zVYdL4B2uNPpNNq7nf1mq9XoNHZf6viHTdvcOei02odtIHfahwedZh3zMIKWDd17xr5uXv4CRrzlLg)


### BamToVariantsWorkFlow

Two workflows a generic germline workflow and a somatic workflow. The somatic workflow is subdivided into freebayes, MuTect2 (gatk4) and LoFreq.

Germline 

<!--
    The workflow should look like:
        flowchart TD

        %% INPUT
        A[Input: SampleConfig + Reference + TargetIntervals]

        %% PREPARE BAMS
        A --\> B{Scatter samples}
        B --\> C[Link BAM]
        B --\> D[Link BAI]

        %% HAPLOTYPECALLER PER SAMPLE + INTERVAL
        C --\> E{Scatter target intervals}
        D --\> E

        E --\> F["HaplotypeCaller (GVCF mode)"]

        %% MERGE PER SAMPLE
        F --\> G[CombineGVCFs per sample]

        %% MERGE ALL SAMPLES
        G --\> H["CombineGVCFs (all samples)"]

        %% GENOTYPING
        H --\> I[GenotypeGVCFs]

        %% OUTPUT
        I --\> J[Final Outputs:\n- Project GVCF\n- Final VCF]
-->

[![](https://mermaid.ink/img/pako:eNpVUluPojAU_itNk0lms2pAvCAPmyBWZINKkJlkl_rQlarsQEtK2d1Z43-fUoe5PJCcc_huPe0FHnhGoQOPBf97OBMhQbLADLO7OxBsoocEMzcNWNVIB-xIWRXU4-yYn8BXENMjFZQdqKoTIk5UBkxS8YcU9f5VIYpR5MYIzN31TgmBfv8bmF92ByIVENRar75iNtd_vDTM2VML3nejRTcKOsmVG4Xb5EeEPDcMUQwi9e3cdRQiFSPYJCh-dEPMPE1Hb15SBwR5l1CZLm6QVhfpcpliuCJVweVzRT1SFIp37z96S1CqHX3BsMuwRrGPPjhjttQCfurx8lfOaEuqQfV2xs9EFfyVqJbia-ZKWX_i3iv7bkEfnH20aQ8fbHzMVpoZpD5lOrEmdsDtQ6IvL9Cg7-kyZ6QA20aqm6wdjFkfRIL_pgcJWp4e3DCq28MePIk8g44UDe3BkoqStC28YAYAhvJMS4qho8qMiCcMMbsqTkXYT87LjiZ4czpD56jWrbqmyoiki5ycBHmHUJZR4fGGSeiY9mSkRaBzgf-gY9mzwcy0xqZhTe3h2LR68FmhLHNgTQ1jYpjDkWWMTfvag_-1rTGYju3J2JoNDWs0MUZDuwdplksu1rdHrt_69QU6eeoF?type=png)](https://mermaid.ai/live/edit#pako:eNpVUluPojAU_itNk0lms2pAvCAPmyBWZINKkJlkl_rQlarsQEtK2d1Z43-fUoe5PJCcc_huPe0FHnhGoQOPBf97OBMhQbLADLO7OxBsoocEMzcNWNVIB-xIWRXU4-yYn8BXENMjFZQdqKoTIk5UBkxS8YcU9f5VIYpR5MYIzN31TgmBfv8bmF92ByIVENRar75iNtd_vDTM2VML3nejRTcKOsmVG4Xb5EeEPDcMUQwi9e3cdRQiFSPYJCh-dEPMPE1Hb15SBwR5l1CZLm6QVhfpcpliuCJVweVzRT1SFIp37z96S1CqHX3BsMuwRrGPPjhjttQCfurx8lfOaEuqQfV2xs9EFfyVqJbia-ZKWX_i3iv7bkEfnH20aQ8fbHzMVpoZpD5lOrEmdsDtQ6IvL9Cg7-kyZ6QA20aqm6wdjFkfRIL_pgcJWp4e3DCq28MePIk8g44UDe3BkoqStC28YAYAhvJMS4qho8qMiCcMMbsqTkXYT87LjiZ4czpD56jWrbqmyoiki5ycBHmHUJZR4fGGSeiY9mSkRaBzgf-gY9mzwcy0xqZhTe3h2LR68FmhLHNgTQ1jYpjDkWWMTfvag_-1rTGYju3J2JoNDWs0MUZDuwdplksu1rdHrt_69QU6eeoF)



### Goals:

- [x] is to get generic alignment working.
- [x] basic variant calling (haplotypecallerGvcf).
- [x] somatic variant calling(MuTect2,...).
- [x] Runs on software installed with Easybuild
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

Easybuild the required modules or use future wrapper easybuild module. 


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
