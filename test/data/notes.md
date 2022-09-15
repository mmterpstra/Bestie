#test/data/
This dir contains test data either validation or for input needed to run the pipeline.

#fastq.gz inputs

Creating fastq inputs for analysis

```
(ml BEDTools/2.27.1-foss-2018a wgsim/a12da33-foss-2018a
	bedtools getfasta -fi /data/umcg-mterpstra/apps/data/ftp.broadinstitute.org/bundle/bundle17jan2020/hg38/Homo_sapiens_assembly38.fasta \
		-bed <( echo -e "chr8\t127736594\t127740958\tNMYC") \
		-fo /dev/stdout > test/data/nmyc.fasta
	 wgsim test/data/nmyc.fasta test/data/fastq/reads_R1.fq test/data/fastq/reads_R2.fq \
		&>/data/umcg-mterpstra/umcg-oncogenetics/git/Bestie/test/data/fastq/reads.wgsim.log
	gzip  test/data/fastq/reads_R1.fq test/data/fastq/reads_R2.fq test/data/nmyc.fasta )
```
