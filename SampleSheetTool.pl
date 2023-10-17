use warnings;
use strict;
use Data::Dumper;
use Getopt::Std;
use POSIX;
#use JSON;
my $use =<<"END";
use: perl $0 [merge|addfile|reformat|reformatmin|stats|uniq|batch|help] 
END
scalar


main();

sub main {
	ToolRunner();
}


sub ToolRunner {
	my $tool = shift @ARGV;
	if(not($tool)){	
		warn $use;
		exit(0);
	}
	if($tool eq 'merge'){
		MergeSampleSheets($tool);
	}elsif($tool eq 'addfile'){
		AddFilesSampleSheet($tool);
	}elsif($tool eq 'jsondump'){
		JsonSampleSheet($tool);
	}elsif($tool eq 'reformat'){
		ReformatSampleSheet($tool);
	}elsif($tool eq 'reformatillumina'){
		ReformatSampleSheetIllumina($tool);
	}elsif($tool eq 'reformatmin'){
                ReformatSampleSheetMin($tool);
        }elsif($tool eq 'stats'){
		CollectStatsSampleSheet($tool);
	}elsif($tool eq 'uniq'){
		UniquefySamplesheets($tool);
	}elsif($tool eq 'valid'){
		ValidateSamplesheet($tool);
	}elsif($tool eq 'batch'){
		BatchSamplesheet($tool);
 	}elsif($tool eq 'split'){
		SplitSamplesheetByProject($tool);
	}elsif($tool eq 'arch'){
		ArchiveDataSamplesheet($tool);
 	}elsif($tool eq 'help'){
		Help($tool);
 	}else{
		die "ERROR tool '$tool' is not correctly specified specified. $use";
	}
}

##############################################################
###Tools
sub Help {
	my $tool = shift @_;
	my $use .=<<"END";
use: perl $0 [merge|addfile|reformat|stats|uniq|help] commands
 toolkit to work with molgenis-compute sample configuration files
 - merge	Merge the tables based on header
	use: perl $0 merge samplesheet_1.csv samplesheet_2.csv [samplesheet_n.csv] > samplesheet.merge.csv

 - addfile	Add fastq filenames to gcc spec based on lane/barcode/sample/read info
	use: perl $0 addfile \
		gccsamplesheet.csv \
		fastqdir > samplesheet.filename.csv

 - reformat	Converts gcc spec to my own spec
	use: perl $0 reformat \
		samplesheet.filename.csv > samplesheet.reformat.csv

 - reformatmin	Converts minimal (tsv with fq1,fq2,fq3,sampleName) spec to my own spec
	use: perl $0 reformatmin \
		samplesheet.filename.csv > samplesheet.reformat.csv

 - stats	Collect statistics on samplesheet
	use: perl $0 stats \

 - uniq 	wip! get uniq record by specifying fields and merge prio by 
		specification order -f specifies the order of fields as csv 
	use: perl $0 uniq \
		samplesheet1.csv samplesheet2.csv \
		-f field1[,field2] >uniq.csv

 - validate	wip! Validate samplesheets using all modules default or only specified modules
	use: perl $0 valid \
		[-m module[,moduleN]]
		samplesheet.csv

 - batch	wip! Batches projectnames based on PROJECTBASE and max INT samples in batches 
		also adds a controlsample in here.
		use: perl $0 $tool -p PROJECTBASE -b INT -c CONTROLSAMPLENAME samplesheet.csv

			
 - help 	This!
	use: find it by `history| tail -n 2` and rtfm
	
END
warn $use
}

######################################################################
#main modules
#to be done FqUniq
#
######################################################################

sub MergeSampleSheets {
	my $tool = shift @_;
	my $use .=<<"END";
use: perl $0 $tool samplesheet_1.csv samplesheet_2.csv [samplesheet_n.csv] > samplesheet.merge.csv 
END
	if(not(scalar(@ARGV))){
		warn $use;
		exit(1);
	}
	my $mergedSamplesheet;
	while(scalar(@ARGV)){
		my $samplesheetfile = shift @ARGV;
                warn "[INFO] Processing samplesheetfile='$samplesheetfile'";

		my $samplesheet = ReadSamplesheet($samplesheetfile);
		#add original file
		for my $sample (@$samplesheet){
			if(not($sample -> {'OrignalCsv'})){
				if(substr(0,1) eq '/'){
					$sample -> {'OrignalCsv'}=$samplesheetfile;
				}else{
					my $pwd=`pwd`;
					chomp($pwd);
					$sample -> {'OrignalCsv'}=$pwd."/".$samplesheetfile;#warn $pwd."/".$samplesheetfile. " ";
				}
			}
		}
		warn "[INFO] imported ".scalar(@$samplesheet)." from samplesheetfile='$samplesheetfile'";
		
		#add to merged samplesheet
		push(@$mergedSamplesheet,@$samplesheet);
	}
	
	print SamplesheetAsString($mergedSamplesheet);
	warn "[INFO] Done";
}

######################################################################

sub AddFilesSampleSheet {
	my $tool = shift @_;
	my $use .=<<"END";
use: perl $0 $tool gccsamplesheet.csv fastqdir > samplesheet.filename.csv
END

	my $samplesheetfile = shift @ARGV;
	my $prefix = shift @ARGV;
	warn "## ".localtime(time())." ## INFO ## $0 init with samplesheetfile='$samplesheetfile';prefix='$prefix'.\n";
	my $targetprefix = shift @ARGV;
	my $samplesheet = ReadSamplesheet($samplesheetfile);
	$samplesheet = AddReadFileNames($samplesheet,$prefix);
	#die Dumper($samplesheet);
	print SamplesheetAsString($samplesheet);
}

##############################################################################

sub ReformatSampleSheet { 
	my $tool = shift @_;
	my $use .=<<"END";
use: perl $0 $tool gccsamplesheetwithfilenames.csv > samplesheet.reformat.csv
END

	my $samplesheetfile = shift @ARGV;
	#my $prefix = shift @ARGV;
	warn "## ".localtime(time())." ## INFO ## $0 init with samplesheetfile='$samplesheetfile'.\n";
	my $samplesheet = ReadSamplesheet($samplesheetfile);
	
	my $reformatTemplate=ReformatHashGccToRnaSeq();
	$samplesheet=ConvertSampleSheet($samplesheet,$reformatTemplate);
	
	#die Dumper($samplesheet);
	print SamplesheetAsString($samplesheet);
}

##############################################################################

sub ReformatSampleSheetMin { 
	my $tool = shift @_;
	my $use .=<<"END";
use: perl $0 $tool minsamplesheetwithfilenames.csv > samplesheet.reformat.csv

note the header should contain fq1/fq2/fq3/fq4 etc
example for creating from scratch:
ls *_R1.fastq.gz| perl -wne 'BEGIN{print "fq1,fq2,sampleName\n"};chomp;print \$_;s/_R1\./_R2./g; print ",\$_"; s/.*-(\\d\\d\\d)_.*/\$1/g; print ",\$_" ;print "\n"'
END

	my $samplesheetfile = shift @ARGV;
	#my $prefix = shift @ARGV;
	warn "## ".localtime(time())." ## INFO ## $0 init with samplesheetfile='$samplesheetfile'.\n";
	my $samplesheet = ReadSamplesheet($samplesheetfile);
	#die Dumper($samplesheet)." ";
	my $samplesheetann = AnnotateSamplesheet($samplesheet);
	#die Dumper($samplesheetann)." ";
	#my $reformatTemplate=ReformatHashGccToRnaSeq();
	#$samplesheet=ConvertSampleSheet($samplesheet,$reformatTemplate);
	
	#die Dumper($samplesheet);
	print SamplesheetAsString($samplesheetann);
}
##############################################################################
sub ReformatSampleSheetIllumina { 
	my $tool = shift @_;
	my $use .=<<"END";
use: perl $0 $tool gccsamplesheetwithfilenames.csv > samplesheet.reformat.csv
END

	my $samplesheetfile = shift @ARGV;
	#my $prefix = shift @ARGV;
	warn "## ".localtime(time())." ## INFO ## $0 init with samplesheetfile='$samplesheetfile'.\n";
	my $samplesheet = ReadSamplesheetIllumina($samplesheetfile);
	
	my $reformatTemplate=ReformatHashIlluminaToSampleSheet();
	$samplesheet=ConvertSampleSheet($samplesheet,$reformatTemplate);
	
	#die Dumper($samplesheet);
	print SamplesheetAsString($samplesheet);
}
#############################################################################
sub JsonSampleSheet { 
	my $tool = shift @_;
	my $use .=<<"END";
Somewhat stupid json dump... wip...
use: perl $0 $tool gccsamplesheetwithfilenames.csv > samplesheet.json
END

	my $samplesheetfile = shift @ARGV;
	#my $prefix = shift @ARGV;
	warn "## ".localtime(time())." ## INFO ## $0 init with samplesheetfile='$samplesheetfile'.\n";
	my $samplesheet = ReadSamplesheet($samplesheetfile);
	#my $samplesheetann = AnnotateSamplesheet($samplesheet);
	#my $reformatTemplate=ReformatHashIlluminaToSampleSheet();
	#$samplesheet=ConvertSampleSheet($samplesheet,$reformatTemplate);
	
	#die Dumper($samplesheet);
	print SamplesheetAsJSONString($samplesheet);
}

######################################################################

sub CollectStatsSampleSheet { 
	my $tool = shift @_;
	my $use .=<<"END";
use: perl $0 $tool gccsamplesheetwithfilenames.csv > samplesheet.reformat.csv
END
	
        my $samplesheetfile = shift @ARGV;
        #my $prefix = shift @ARGV;
        warn "## ".localtime(time())." ## INFO ## $0 init with samplesheetfile='$samplesheetfile'.\n";
        my $samplesheet = ReadSamplesheet($samplesheetfile);

        my $stats=StatsCollecterSampleSheet($samplesheet);


	exit();
}

########################################################################

sub UniquefySamplesheets { 
	my $tool = shift @_;
	my $use .=<<"END";
use: perl $0 $tool -f field, samplesheet.csv samplesheet2.csv > samplesheet.uniq.csv
END
	my $opts;%{$opts}=();
	getopts('f:', $opts);
	my $mergedSamplesheet;
	
	if(not(scalar(@ARGV))){
		warn $use;
		exit(1);
	}

	while(scalar(@ARGV)){
		my $samplesheetfile = shift @ARGV;
                warn "[INFO] Processing samplesheetfile='$samplesheetfile'";

		my $samplesheet = ReadSamplesheet($samplesheetfile);
		#add original file
		for my $sample (@$samplesheet){
			if(not($sample -> {'OrignalCsv'})){
				if(substr(0,1) eq '/'){
					$sample -> {'OrignalCsv'}=$samplesheetfile;
				}else{
					my $pwd=`pwd`;
					chomp($pwd);
					$sample -> {'OrignalCsv'}=$pwd."/".$samplesheetfile;#warn $pwd."/".$samplesheetfile. " ";
				}
			}
		}
		warn "[INFO] imported ".scalar(@$samplesheet)." from samplesheetfile='$samplesheetfile'";
		
		#add to merged samplesheet
		push(@$mergedSamplesheet,@$samplesheet);
	}
	
	#subroutines
	my @fields;
	my $fields;%{$fields} =();

	my $samplesheetNew;
	@$samplesheetNew=();
	warn "[INFO] start uniquefy\n";
	
	if($opts -> {'f'}){
		my $h;%{$h}=();
		
		#collect samples by fields
		@fields = split(',',$opts -> {'f'});
		for my $sample (@$mergedSamplesheet){
			my $ref = $h;
			for my $field (@fields){
				my $val = $sample -> {$field};
				if(not(defined($ref -> {$val}))){
					%{$ref -> {$val}} = ();
				}
				#move the hash ref a step down the construct
				$ref = $ref -> {$val};
			}
			push(@{$ref -> {'samples'}},$sample);
		}
		#Iterate $h leaves which are arrays of samples and prio samples
		#https://stackoverflow.com/questions/2363142/how-to-iterate-through-hash-of-hashes-in-perl
		my $samplelist;
		@{$samplelist}=();
		HashWalkBestSamples($h, [],$samplelist, \&CollectBestSamples);		
		
		warn "[INFO] end uniquefy\n";
		print SamplesheetAsString($samplelist);
		#for my $mod (@mods){
		#	die "Invalid mod $mod this should be 'controls', 'fastq', chars, cols and/or 'id'" if(not($mod eq 'controls' || $mod eq 'fastq' || $mod eq 'id' || $mod eq 'chars' || $mod eq 'cols'));
		#	$mods -> {$mod} ++ if($mod eq 'chars' || $mod eq 'cols');
		#}
	}else{
		warn "Specify field(s) -f field1[,field2]\n";
		warn $use;
	
		exit(1);		
		#run all mods
		#$mods -> {'cols'} ++;
		#$mods -> {'chars'} ++;
	}
	
        my $samplesheetfile = shift @ARGV;
        #my $prefix = shift @ARGV;
        warn "## ".localtime(time())." ## INFO ## $0 end.\n";


	exit();
}

##########################################################################

sub ValidateSamplesheet {
	my $tool = shift @_;
	my $opts;%{$opts}=();
	getopts('m:', $opts);
	my $mergedSamplesheet;
	my $use .=<<"END";
WIP use: perl $0 $tool [-m MODULE,[MODULE]] samplesheet.csv
check for valid samplesheet files using different features with modules. defaults to all modules and possible vals:'controls', 'fastq', chars, cols and/or 'id.
END

	#subroutines
	my @subs = ("controls","fastq","id","cols");
	my $mods;%{$mods} =();

	my @mods;
	if($opts -> {m}){
		@mods = split(',',$opts -> {m});
		for my $mod (@mods){
			die "Invalid mod $mod this should be 'controls', 'fastq', chars, cols and/or 'id'" if(not($mod eq 'controls' || $mod eq 'fastq' || $mod eq 'id' || $mod eq 'chars' || $mod eq 'cols'));
			$mods -> {$mod} ++ if($mod eq 'chars' || $mod eq 'cols');
		}
	}else{
		#run all mods
		$mods -> {'cols'} ++;
		$mods -> {'chars'} ++;
	}
	
	while(scalar(@ARGV)){
		my $samplesheetfile = shift @ARGV;
		
		my $samplesheet = ReadSamplesheetWithValidation($samplesheetfile,$mods);
		push(@$mergedSamplesheet,@$samplesheet);
	}
	if($opts -> {m}){
		for my $mod (@mods){
			
			#ValidateHeader# eq 'head' some of these are in teh old file incoded in the reader and should be separated to toggle
			#ValidateColvalues#eq 'cols'
			#ValidateSpecialChars#eq 'chars'
			ValidateControlSampleNames($mergedSamplesheet) if($mod eq 'controls');
			ValidateFqFiles($mergedSamplesheet) if($mod eq 'fastq');
			ValidateInternalId($mergedSamplesheet) if($mod eq 'id');
			ValidateSequence($mergedSamplesheet) if($mod eq 'seq');
		}
	}else{
		ValidateControlSampleNames($mergedSamplesheet);
		ValidateFqFiles($mergedSamplesheet);
		ValidateInternalId($mergedSamplesheet);
		ValidateSequence($mergedSamplesheet);#ILLUMINA,SLX,SOLEXA,SOLID,454,LS454,COMPLETE,PACBIO,IONTORRENT,CAPILLARY,HELICOS,UNKNOWN
	}
	print SamplesheetAsString($mergedSamplesheet);
}

##########################################################################

sub BatchSamplesheet {
	my $tool = shift @_;
	my $opts;%{$opts}=();
	getopts('p:b:c:m:', $opts);
	my $samplesheetfile = shift @ARGV;
	my $use .=<<"END";
WIP use: perl $0 $tool -p PROJECTBASE -b INT -c CONTROLSAMPLENAME -m MOARCONTROLS samplesheet.csv
Batches projectnames based on PROJECTBASE and about INT samples in batches also adds a CONTROLSAMPLE in here. When MOARCONTROLS is a comma separated list of samplenames these get added as additional controlsamples.
END
	die "ERROR: At least one of these missing -p PROJECTBASE -b INT -c CONTROLSAMPLENAME -m MOARCONTROLS. $use." if(not( $opts -> {'p'} && $opts -> {'b'} && $opts -> {'c'}&& $opts -> {'m'} ));	
	
	my $samplesheet = ReadSamplesheet($samplesheetfile);
	my $batchedss = Batcher($samplesheet,$opts -> {'p'},$opts -> {'b'},$opts -> {'c'},$opts -> {'m'});
	
	print SamplesheetAsString($batchedss);
}

################################

sub SplitSamplesheetByProject {
	my $tool = shift @_;
	my $opts;%{$opts}=();
	getopts('p:b:c:m:', $opts);
	my $samplesheetfile = shift @ARGV;
	my $use .=<<"END";
WIP use: perl $0 $tool -p PREFIX samplesheet.csv
Splits samplesheets based on projects and will put them in prefix_projectname.csv
END
	die "ERROR: At least one of these missing -p PREFIX . $use." if(not( $opts -> {'p'}));	
	
	my $samplesheet = ReadSamplesheet($samplesheetfile);
	my @splitss = SplitByProject($samplesheet);
	
	for my $ss (@splitss){
		warn "INFO: Creating ".$opts -> {'p'}."_".$ss -> [0] -> {'project'}.".csv ...";
		open(my $fh, ">", $opts -> {'p'}."_".$ss -> [0] -> {'project'}.".csv"	);
		print $fh SamplesheetAsString($ss);
		close $fh;
	}
}

################################

sub ArchiveDataSamplesheet {
	my $tool = shift @_;
	my $opts;%{$opts}=();
	getopts('i:a:', $opts);
	my $samplesheetfile = shift @ARGV;
	my $use .=<<"END";
WIP use: perl $0 $tool -i inputdir -a archivedir samplesheet.csv
	Gets all samples/fastqs using the samplesheet and puts em into the archive 
	dir with the samplesheets. Also does adding samplesheets , checking md5sums and making 
	the struture read only. 
END
	die "ERROR: At least one of these missing -p PREFIX . $use." if(not( $opts -> {'p'}));	
	
	my $samplesheet = ReadSamplesheet($samplesheetfile);
	#my @splitss = SplitByProject($samplesheet);
	
	for my $s (@{$samplesheet}){
		warn "INFO: Checking and archiving ".$s -> {sampleName}." to ".$opts -> {'a'} . "\n";
		open(my $fh, ">", $opts -> {'a'}."/".".samplesheet.csv"	);
		print $fh SamplesheetAsString($s);
		close $fh;
	}
}


#################################
#Subroutines
####
##
#
#

sub ReadSamplesheet {
	my $samplesheetf = shift @_;
	my $samplesheet;
	open( my $samplesheeth,'<', $samplesheetf) or die "Cannot read samplesheet file:'$samplesheetf'";
	$_=<$samplesheeth>;
	chomp;
	my @h = CommaseparatedSplit($_);
	#die Dumper(\@h);
	while(<$samplesheeth>){
		chomp;
		my @c = CommaseparatedSplit($_);
		#die Dumper(\@c);
		if(scalar(@c)==scalar(@h)){
			my %d;
			my $i=0;
			map{$d{$_}=$c[$i]; $i++}(@h);
			$c[$i]=join(",",@h);
			#ReadFileNameConstructor(\%d);
			push(@$samplesheet,\%d);
		}else{
			die "## ERROR ## Header is not of equal length compared sample line header:".Dumper(\@h)."columns".Dumper(\@c)." ";
		}
	}
	return $samplesheet;
}
sub AnnotateSamplesheet {
	my $samplesheet = shift @_;
	my $reformatTemplate = shift @_;
	
	#die Dumper($samplesheet)."  here";
	my $newSamplesheet;
	my $id = 1;
	for my $sample (@$samplesheet){
		my $newSample;
		$newSample -> {'sampleName'} = $sample -> {'sampleName'};
		#missing:internalId,samplePrep,sequencer,sequencerId,barcode,sequencingStartDate,project,
		#	controlSampleName,reads3FqGz
		$newSample -> {'reads1FqGz'}="something went wrong using samplesheet tool.";
		$newSample -> {'internalId'}="$id";$id++;
		$newSample -> {'samplePrep'}="oneLibraryPrepPerSampleAssumed";
		$newSample -> {'sequencer'}="illumina";#assumed
		my $date;@{$date} = CmdRunner('date --iso-8601=date');
		chomp($date -> [-1]);
		$newSample -> {'sequencingStartDate'}=$date -> [0];
		$newSample -> {'project'}="projectNameHere";#lazy find and replace for this keyword in generatescrpts
		$newSample -> {'controlSampleName'}=$sample -> {'sampleName'};
		$newSample -> {'reads2FqGz'}="";#always present but empty when not found
		$newSample -> {'reads3FqGz'}="";
		$newSample -> {'barcode'} = "NNNNNN"; 
		my @files;
		map{push(@files,$sample -> {$_})if($_ =~ m/^fq\d/)}(keys(%{$sample}));
		#die Dumper(@files);
		#i put in a sort but this might still give problems based on the sorting of filenames....
		for my $file (sort(@files)){
			chomp $file;
			#Add defaults here when not present in fastq files
			#unreadble code.I hope it works!
			if($file =~ m/.txt.gz/){#This is very experimental should work for single end sequencing
				$newSample -> {'reads1FqGz'} = $file;
				my $fqhead;@{$fqhead} = CmdRunner('zcat '.$newSample -> {reads1FqGz}.'|head -n 1');
				my @fqheadsplit= split /[\/_\-: \@#]/,($fqhead -> [0]);
				chomp($fqheadsplit[-1]);
				#@HWI-ST001_0001:1:1111:112345:12345#NNNNNN/1 for read end 1
				$newSample -> {'run'} = $fqheadsplit[3];
				$newSample -> {'sequencerId'} = $fqheadsplit[2];
				$newSample -> {'flowcellId'} = "nofcid";
				$newSample -> {'seqType'} = "illumina";
				$newSample -> {'lane'} = $fqheadsplit[4];#this should create an uniq list of lanes
				$newSample -> {'barcode'}=$fqheadsplit[-2] if($fqheadsplit[-2] =~ m/[ATCGN\+]{6,}/);
				$newSample -> {'reads1FqGzMd5'}=Md5Sum($newSample -> {"reads1FqGz"});
			}if($file =~ m/_R1\./){
				$newSample -> {'reads1FqGz'} = $file;
				#zcat $file |head -n 1 should contain @M00000:999:000000000-FLOWCELLIDD:TILE:12345:23456:4567 1:N:0:GTGATTCC+TATAGCCT
				my $fqhead;@{$fqhead} = CmdRunner('zcat '.$newSample -> {reads1FqGz}.'|head -n 1');
				my @fqheadsplit= split /[: \@]/,($fqhead -> [0]);
				chomp($fqheadsplit[-1]);
				$newSample -> {'run'} = $fqheadsplit[2];
				$newSample -> {'sequencerId'} = $fqheadsplit[1];
				$newSample -> {'flowcellId'} = $fqheadsplit[3];
				$newSample -> {'seqType'} = $fqheadsplit[3];
				$newSample -> {'lane'} = $fqheadsplit[4];#this should create an uniq list of lanes
				$newSample -> {'barcode'}=$fqheadsplit[-1] if($fqheadsplit[-1] =~ m/[ATCGN\+]{6,}/);
				$newSample -> {'reads1FqGzMd5'}=Md5Sum($newSample -> {"reads1FqGz"});
			}elsif($file =~ m/_R2\./){
				$newSample -> {'reads2FqGz'} = $file;
				my $fqhead;@{$fqhead} = CmdRunner('zcat '.$newSample -> {reads2FqGz}.'|head -n 1');
				my @fqheadsplit= split /[: \@]/,($fqhead -> [0]);
				$newSample -> {'run'} = $fqheadsplit[2];
				$newSample -> {'flowcellId'} = $fqheadsplit[3];
				$newSample -> {'seqType'} = $fqheadsplit[3];
				$newSample -> {'lane'} = $fqheadsplit[4];
                                $newSample -> {'reads2FqGzMd5'}=Md5Sum($newSample -> {"reads2FqGz"});

			}elsif($file =~ m/_R3\./){#nugene/umi probs
				$newSample -> {'reads3FqGz'} = $newSample -> {reads2FqGz};	
				$newSample -> {'reads3FqGzMd5'} = $newSample -> {reads2FqGzMd5};	
				$newSample -> {'reads2FqGz'} = $file;
				my $fqhead;@{$fqhead} = CmdRunner('zcat '.$newSample -> {reads2FqGz}.'|head -n 1');
				my @fqheadsplit= split /[: \@]/,($fqhead -> [0]);
				$newSample -> {'run'} = $fqheadsplit[2];
				$newSample -> {'flowcellId'} = $fqheadsplit[3];
				$newSample -> {'seqType'} = $fqheadsplit[3];
				$newSample -> {'lane'} = $fqheadsplit[4];#this should create an uniq list of lanes
                                $newSample -> {'reads2FqGzMd5'}=Md5Sum($newSample -> {"reads2FqGz"});
			
			}elsif($file =~ m/_umi\./){#umi probes literal no annoying R2 as umi and R2=R3
				$newSample -> {'reads3FqGz'} = $file;
				my $fqhead;@{$fqhead} = CmdRunner('zcat '.$newSample -> {reads2FqGz}.'|head -n 1');
				my @fqheadsplit= split /[: \@]/,($fqhead -> [0]);
				$newSample -> {'run'} = $fqheadsplit[2];
				$newSample -> {'flowcellId'} = $fqheadsplit[3];
				$newSample -> {'seqType'} = $fqheadsplit[3];
				$newSample -> {'lane'} = $fqheadsplit[4];#this should create an uniq list of lanes
                                $newSample -> {'reads3FqGzMd5'}=Md5Sum($newSample -> {"reads3FqGz"});
			
			}elsif($file =~ m/_R4\./){
				#Future proofing
				$newSample -> {'reads4FqGz'} = $file;
				$newSample -> {'reads4FqGzMd5'} = Md5Sum($newSample -> {"reads4FqGz"});
				my $fqhead;@{$fqhead} = CmdRunner('zcat '.$newSample -> {reads4FqGz}.'|head -n 1');
				my @fqheadsplit= split /[: \@]/,($fqhead -> [0]);
				$newSample -> {'run'} = $fqheadsplit[2];
				$newSample -> {'flowcellId'} = $fqheadsplit[3];
				$newSample -> {'seqType'} = $fqheadsplit[3];
				$newSample -> {'lane'} = $fqheadsplit[4];#this should create an uniq list of lanes
			}
		}
		#die Dumper($newSample);
		push (@$newSamplesheet, $newSample);
	}
	return $newSamplesheet;
}
sub ReadSamplesheetIllumina {
	my $samplesheetf = shift @_;
	my $samplesheet;
	open( my $samplesheeth,'<', $samplesheetf) or die "Cannot read samplesheet file:'$samplesheetf'";
	my $headerandsettings;
	my $line = " ";
	while(defined($line = <$samplesheeth>) and not( $line =~ m/\[Data\]/)){
		#warn "line " . $line; 
		$line =~  s/\r?\n\z//;
		#$line =~  s/\r|\n|\z//;		
		chomp $line;
		next if($line =~ m/^\[/);
		if($line =~ m/\,/){
			my @c = CommaseparatedSplit($line);
			#warn "###Commaline     " . $line ." \n". Dumper(\@c). " ";
			#$c[0] =~  s/\r|\n|\z//;		
			$headerandsettings -> {$c[0]} = $c[1];
		}
	}
	#die "Headerandsettings\n".Dumper($headerandsettings)." ";
	$line = <$samplesheeth>;
	$line =~  s/\r?\n\z//;
	chomp $line;
	my @h = CommaseparatedSplit($line);
		
	#die Dumper(\@h);
	while(<$samplesheeth>){
		s/\r?\n\z//;		
		chomp;
		my @c = CommaseparatedSplit($_);
		#die Dumper(\@c);
		if(scalar(@c)==scalar(@h)){
			my %d = %{$headerandsettings};
			$d{'sequencer'}= 'illumina';
			$d{'line'}=$.;
			my $i=0;
			map{$d{$_}=$c[$i]; $i++}(@h);
			$c[$i]=join(",",@h);
			die "ERROR Invalid sample id: '".$d{'Sample_ID'}."'. This should not contain spaces." if($d{'Sample_ID'} =~ m/ /);
			#ReadFileNameConstructor(\%d);
			my $ret;
			@{$ret} = CmdRunner('ls "$(dirname '."'$samplesheetf'".')"/'.$d{'Sample_Project'}.'/'.$d{'Sample_ID'}.'_*_R[1234]_*.fastq.gz');
			if(scalar(@{$ret})){
				my @lanes;
				for my $file (@{$ret}){
					chomp $file;
					#unreadble code.I hope it works!
					if($file =~ m/_R1_0/){
						$d{'reads1FqGz'} = $file;
						#zcat $file |head -n 1 should contain @M00000:999:000000000-FLOWCELLIDD:TILE:12345:23456:4567 1:N:0:GTGATTCC+TATAGCCT
						my $fqhead;@{$fqhead} = CmdRunner('zcat '.$d{reads1FqGz}.'|head -n 1');
						my @fqheadsplit= split /[: \@]/,($fqhead -> [0]);
						$d{'run'} = $fqheadsplit[2];
						$d{'flowcellId'} = $fqheadsplit[3];
						$d{'seqType'} = $fqheadsplit[3];
						push @lanes,  $fqheadsplit[4];#this should create an uniq list of lanes
					}elsif($file =~ m/_R2_0/){
						$d{'reads2FqGz'} = $file;
						my $fqhead;@{$fqhead} = CmdRunner('zcat '.$d{reads1FqGz}.'|head -n 1');
						my @fqheadsplit= split /[: \@]/,($fqhead -> [0]);
						$d{'run'} = $fqheadsplit[2];
						$d{'flowcellId'} = $fqheadsplit[3];
						$d{'seqType'} = $fqheadsplit[3];
					}elsif($file =~ m/_R3_0/){#nugene/umi probs
						$d{'reads3FqGz'} = $d{reads2FqGz};		
						$d{'reads2FqGz'} = $file;
						my $fqhead;@{$fqhead} = CmdRunner('zcat '.$d{reads1FqGz}.'|head -n 1');
						my @fqheadsplit= split /[: \@]/,($fqhead -> [0]);
						$d{'run'} = $fqheadsplit[2];
						$d{'flowcellId'} = $fqheadsplit[3];
						$d{'seqType'} = $fqheadsplit[3];
					}elsif($file =~ m/_R4_0/){
						#Future proofing
						$d{'reads4FqGz'} = $file;
						my $fqhead;@{$fqhead} = CmdRunner('zcat '.$d{reads1FqGz}.'|head -n 1');
						my @fqheadsplit= split /[: \@]/,($fqhead -> [0]);
						$d{'run'} = $fqheadsplit[2];
						$d{'flowcellId'} = $fqheadsplit[3];
						$d{'seqType'} = $fqheadsplit[3];
					}
				}
				for my $lane (@lanes){
					$d{'lane'} = $lane;
					for my $fq ('reads1FqGz','reads2FqGz','reads3FqGz'){
						my $prettylane = '_L'.sprintf("%03d", $lane).'_';
						if(defined($d{$fq})){
							$d{$fq} =~ s/_L\d\d\d_/$prettylane/ge;
							$d{$fq."Md5"} = Md5Sum($d{$fq});
						}if(not(defined($d{$fq}))){
							$d{$fq} = '';
						} 
					}
					
					push(@$samplesheet,\%d);
				}
			}else{ 
				push(@$samplesheet,\%d);
			}
		}else{
			die "## ERROR ## Header is not of equal length compared sample line header:".Dumper(\@h)."columns".Dumper(\@c)." ";
		}
	}
	#die Dumper($samplesheet)." ";
	return $samplesheet;
}
sub CommaseparatedSplit {
	my $string=pop @_;
	#needs to be fixed for citation marks!
	warn "Line contains citation marks: this is currently not supported!!. I hope this works. Line=$_" if($string =~ /"|'/);
	my $i = index($string,",");
	if( $i > -1){
		push(@_,substr($string,0,$i));
		push(@_,substr($string,$i+1));
		@_ = CommaseparatedSplit( @_ );
	}else{
		push(@_,$string);
		return @_;
	}
}

sub SamplesheetAsString {
	my $self = shift @_;
	my $string = '';
	#get header values;
	my %h;
	for my $sample (@$self){
		for my $key (keys(%$sample)){
			$h{$key}++;
		}
	}
	my @h = sort {$b cmp $a} (keys(%h));
	$string.=join(",",@h)."\n";
	#warn scalar(@$self);
	for my $sample (@$self){
		my @c;
		for my $h (@h){
			if($$sample{$h}){
				push (@c,$$sample{$h});
			}else{
				push (@c,"");
			}
		}
		$string.=join(",",@c)."\n";
	}
	return $string;
}

sub SamplesheetAsJSONString {
	my $self = shift @_;
	my $string = '';
	#get header values;
	my %h;
	for my $sample (@$self){
		for my $key (keys(%$sample)){
			$h{"KEYS"} -> {$key}++ if not $key eq "sampleName";
			push(@{$h{"SAMPLENAMES"} -> {$sample -> {'sampleName'}}},$sample) if $key eq "sampleName";
		}
	}
	
	$string.= '{"samples": [';
	
	#my @h = sort {$b cmp $a} (keys(%h));
	#$string.=join(",",@h)."\n";
	#warn scalar(@$self);
	my $idx = 0;
	for my $samplename (keys(%{$h{"SAMPLENAMES"}})){
		#die Dumper(%h);

		$string.= '    {"name": "'.$h{"SAMPLENAMES"}->{$samplename} ->  [0] -> {'sampleName'}.'",
			"threeLetterName": "MS1",
			"readgroups" : ['."\n";
		my $rgidx = 0;
		for my $rg (@{$h{"SAMPLENAMES"} -> { $samplename }}){
		#	my @c;
		#	for my $h (@h){
		#		if($$sample{$h}){
		#			push (@c,$$sample{$h});
		#		}else{
		#			push (@c,"");
		#		}
		#	}
		#	$string.=join(",",@c)."\n";

			$string.= '    				{"identifier": "'.$rg -> {'sampleName'}.'_'.
						$rg -> {'lane'}.'_'.$rg -> {'barcode'}.'",
					"lane": "'.$rg -> {'lane'}.'",
					"barcode2": "NNNNNN",
					"fastq1": "'.$rg -> {'reads1FqGz'}.'",'."\n";
			$string .= '                "fastq2": "'.$rg -> {'reads2FqGz'}.'",'."\n" if($rg -> {'reads2FqGz'} ne "");
					$string .= '					"date": "'.$rg -> {'sequencingStartDate'}.'",
					"barcode1": "'.$rg -> {'barcode'}.'",
					"run": "'.$rg -> {'run'}.'",
					"platform_model" : "HiSeq2xxx",
					"sequencer":"'.$rg -> {'sequencerId'}.'",
					"platform": "'.$rg -> {'sequencer'}.'",
					"flowcell":"'.$rg -> {'flowcellId'}.'",
					"sequencing_center": "undef"
				}';
			$string = $string .',' if(($rgidx+1) < scalar(@{$h{"SAMPLENAMES"} -> { $samplename }}));
			$string .= "\n";
			$rgidx++;
		}

		$string .= '       ]'."\n".'    }';
		$string = $string .',' if(($idx+1) < scalar(keys(%{$h{"SAMPLENAMES"}})));
		$string .= "\n";
		$idx++;
	}
	$string.= ']}'."\n";
	
	return $string;
}

sub ReadBaseNameBuilder {
	my $self = shift @_;
	
	#140430_SN163_0549_AH78TFADXX_L1_GAGTAGAG_2
	#die "columns PREFIX and TARGETPRFIX need to be present. specify on commandline!" if($$self{'PREFIX'} && $$self{'TARGETPRFIX'});
	die "Column sequencingStartDate needs to be present. Specify in csv file!".Dumper($self) if(not(exists($$self{'sequencingStartDate'})));
	die "Column sequencer needs to be present. Specify in csv file!".Dumper($self) if(not(exists($$self{'sequencer'})));
	die "Column run to be present. Specify in csv file!".Dumper($self) if(not(exists($$self{'run'})));
	die "Column flowcell needs to be present. Specify in csv file!".Dumper($self) if(not(exists($$self{'flowcell'})));
	die "Column lane needs to be present. Specify in csv file!".Dumper($self) if(not(exists($$self{'lane'})));
	die "Column barcode needs to be present. Specify in csv file!".Dumper($self) if(not(exists($$self{'barcode'})) && not(exists($$self{'barcode1'})));
	
	
	my @fqBases=();
	
	#GCC version of fastq files
	if(exists($$self{'barcode'})){
		push(@fqBases,$$self{'sequencingStartDate'}."_".$$self{'sequencer'}."_".$$self{'run'}."_".$$self{'flowcell'}."_L".$$self{'lane'}."_".$$self{'barcode'});
		push(@fqBases,$$self{'sequencingStartDate'}."_".$$self{'sequencer'}."_".'0'x(4- length($$self{'run'})).$$self{'run'}."_".$$self{'flowcell'}."_L".$$self{'lane'}."_".$$self{'barcode'});
	}elsif(exists($$self{'barcode1'}) && exists($$self{'barcode2'}) && $$self{'barcode1'} ne "NA" &&  $$self{'barcode2'} ne "NA"){
		push(@fqBases,$$self{'sequencingStartDate'}."_".$$self{'sequencer'}."_".$$self{'run'}."_".$$self{'flowcell'}."_L".$$self{'lane'}."_".$$self{'barcode1'}."-".$$self{'barcode2'});
		push(@fqBases,$$self{'sequencingStartDate'}."_".$$self{'sequencer'}."_".'0'x(4- length($$self{'run'})).$$self{'run'}."_".$$self{'flowcell'}."_L".$$self{'lane'}."_".$$self{'barcode1'}."-".$$self{'barcode2'});
	}elsif(exists($$self{'barcode1'}) &&  $$self{'barcode1'} ne "NA"){
		push(@fqBases,$$self{'sequencingStartDate'}."_".$$self{'sequencer'}."_".$$self{'run'}."_".$$self{'flowcell'}."_L".$$self{'lane'}."_".$$self{'barcode1'});
		push(@fqBases,$$self{'sequencingStartDate'}."_".$$self{'sequencer'}."_".'0'x(4- length($$self{'run'})).$$self{'run'}."_".$$self{'flowcell'}."_L".$$self{'lane'}."_".$$self{'barcode1'});
	}else{
		push(@fqBases,$$self{'sequencingStartDate'}."_".$$self{'sequencer'}."_".$$self{'run'}."_".$$self{'flowcell'}."_L".$$self{'lane'});
		push(@fqBases,$$self{'sequencingStartDate'}."_".$$self{'sequencer'}."_".'0'x(4- length($$self{'run'})).$$self{'run'}."_".$$self{'flowcell'}."_L".$$self{'lane'});
	}

	#my version
	push(@fqBases,"lane".$$self{'lane'}."_".$$self{'barcode'}.'_S*_L00'.$$self{'lane'}."_R1_001");
	
	#default illumina version
	#samplename_S5_L001_I1_001.fastq.gz
	push(@fqBases,$$self{'externalSampleID'}.'_S*_L00'.$$self{'lane'}."_R1_001");
	
	warn "looking for ". join(', ', @fqBases)." \n";

	#lane4_TGTTCCGT_S128_L004_R2_001.fastq.gz
	return @fqBases;
	
	
	#return $self;
}

sub AddReadFileNames {
	my $samplesheet = shift @_;
	my $prefix = shift @_;
	#my $targetprefix = shift @_;
	
	my $samplesheetNew;
	@{$samplesheetNew}=();
	
	my $scount=0;
	while($scount < scalar(@{$samplesheet})){
		my $sample = $samplesheet->[$scount];
		my $samplePushed=0;
		my @fqbases = ReadBaseNameBuilder($sample);
		warn scalar @fqbases;
		for my $fqbase (@fqbases){
			
			#warn $prefix."/".$fqbase."";
			#warn $prefix."/".$fqbase."_1.fq.gz" if(-e $prefix."/".$fqbase."_1.fq.gz");
			my $fqbase2 = $fqbase;
			my $fqbase3 = $fqbase;
			my $fqbaseI1 = $fqbase;
			if(index($fqbase,"_R1_")){
				substr($fqbase2,index($fqbase2,"_R1_"),4,"_R2_");
				substr($fqbase3,index($fqbase3,"_R1_"),4,"_R3_");
				substr($fqbaseI1,index($fqbaseI1,"_R1_"),4,"_I1_");
			}
			my @exts=('.fq.gz', '.fastq.gz');
			#glob or globble?
			for my $ext (@exts){
				#warn "looking for $fqbase\n$fqbase2\n$fqbase3\n";
				if( -e $prefix."/".$fqbase.$ext){
					
					$sample -> {"reads1FqGz"}=$prefix."/".$fqbase.$ext;
					$sample -> {"reads1FqGzMd5"}=Md5Sum($sample -> {"reads1FqGz"});
					warn "## ".localtime(time())." ## INFO ## Single end read detected. Stored as read1 = ".$$sample{'reads1FqGz'}. "; md5sum = ".$sample -> {'reads1FqGzMd5'};
					
				}elsif((-e $prefix."/".$fqbase."_1".$ext && ! -e $prefix."/".$fqbase."_2".$ext )){
					
					$sample -> {"reads1FqGz"}=$prefix."/".$fqbase."_1".$ext;
					$sample -> {"reads1FqGzMd5"}=Md5Sum($sample -> {"reads1FqGz"});
					warn "## ".localtime(time())." ## WARN ## Paired end read file name detected but only one end found thus single end assumed. Stored as read1 = ".$$sample{'reads1FqGz'};
					
				}elsif(( ! -e $prefix."/".$fqbase."_1".$ext && -e $prefix."/".$fqbase."_2".$ext )){
					
					$sample -> {"reads1FqGz"}=$prefix."/".$fqbase."_2".$ext;
					$sample -> {"reads1FqGzMd5"}=Md5Sum($sample -> {"reads1FqGz"});
					warn "## ".localtime(time())." ## WARN ## Paired end read file name detected but only one end found thus single end assumed. Stored as read1 = ".$$sample{'reads1FqGz'};
					
				}elsif(( -e $prefix."/".$fqbase."_1".$ext && -e $prefix."/".$fqbase."_2".$ext )){
					
					$sample -> {"reads1FqGz"}=$prefix."/".$fqbase."_1".$ext;
					$sample -> {"reads1FqGzMd5"}=Md5Sum($sample -> {"reads1FqGz"});

					$sample -> {"reads2FqGz"}=$prefix."/".$fqbase."_2".$ext;
					$sample -> {"reads2FqGzMd5"}=Md5Sum($sample -> {"reads2FqGz"});

					warn "## ".localtime(time())." ## INFO ## Paired end read files detected. Stored as read1 = ".$$sample{'reads1FqGz'}." read2 = ".$$sample{'reads2FqGz'};
					#warn Dumper($sample);

				#illumina raw output rules
				#pe with index
				}elsif(index($fqbase,"_R1_") && ( -e glob($prefix."/".$fqbase.$ext) && -e glob($prefix."/".$fqbase2.$ext) && -e glob($prefix."/".$fqbase3.$ext) && -e glob($prefix."/".$fqbaseI1.$ext)) ){
					
					$sample -> {"reads1FqGz"}=glob($prefix."/".$fqbase.$ext);
					$sample -> {"reads1FqGzMd5"}=Md5Sum($sample -> {"reads1FqGz"});

					$sample -> {"reads2FqGz"}=glob($prefix."/".$fqbase3.$ext);
					$sample -> {"reads2FqGzMd5"}=Md5Sum($sample -> {"reads2FqGz"});

					$sample -> {"reads3FqGz"}=glob($prefix."/".$fqbase2.$ext);
					$sample -> {"reads3FqGzMd5"}=Md5Sum($sample -> {"reads3FqGz"});

					$sample -> {"reads1IndexFqGz"}=glob($prefix."/".$fqbaseI1.$ext);
					$sample -> {"reads1IndexFqGzMd5"}=Md5Sum($sample -> {"reads1IndexFqGz"});

					warn "## ".localtime(time())." ## INFO ## Paired end read + index files detected. Stored as read1 = ".$$sample{'reads1FqGz'}." read2 = ".$$sample{'reads2FqGz'};
					#warn Dumper($sample);
				#pe without index
				}elsif(index($fqbase,"_R1_") && ( -e glob($prefix."/".$fqbase.$ext) && -e glob($prefix."/".$fqbase2.$ext) && -e glob($prefix."/".$fqbase3.$ext)) ){
					
					$sample -> {"reads1FqGz"}=glob($prefix."/".$fqbase.$ext);
					$sample -> {"reads1FqGzMd5"}=Md5Sum($sample -> {"reads1FqGz"});

					$sample -> {"reads2FqGz"}=glob($prefix."/".$fqbase3.$ext);
					$sample -> {"reads2FqGzMd5"}=Md5Sum($sample -> {"reads2FqGz"});

					$sample -> {"reads3FqGz"}=glob($prefix."/".$fqbase2.$ext);
					$sample -> {"reads3FqGzMd5"}=Md5Sum($sample -> {"reads3FqGz"});

					$sample -> {"reads1IndexFqGz"}=glob($prefix."/".$fqbaseI1.$ext);
					$sample -> {"reads1IndexFqGzMd5"}=Md5Sum($sample -> {"reads1IndexFqGz"});

					warn "## ".localtime(time())." ## INFO ## Paired end read files detected. Stored as read1 = ".$$sample{'reads1FqGz'}." read2 = ".$$sample{'reads2FqGz'};
					#warn Dumper($sample);
				#se with index
				}elsif(index($fqbase,"_R1_") && ( -e glob($prefix."/".$fqbase.$ext) && -e glob($prefix."/".$fqbase2.$ext) && -e glob($prefix."/".$fqbaseI1.$ext)) ){
					
					$sample -> {"reads1FqGz"}=glob($prefix."/".$fqbase.$ext);
					$sample -> {"reads1FqGzMd5"}=Md5Sum($sample -> {"reads1FqGz"});

					$sample -> {"reads2FqGz"}="";
					$sample -> {"reads3FqGz"}=glob($prefix."/".$fqbase2.$ext);
					$sample -> {"reads3FqGzMd5"}=Md5Sum($sample -> {"reads3FqGz"});

					$sample -> {"reads1IndexFqGz"}=glob($prefix."/".$fqbaseI1.$ext);
					$sample -> {"reads1IndexFqGzMd5"}=Md5Sum($sample -> {"reads1IndexFqGz"});

					warn "## ".localtime(time())." ## INFO ## Single end + index read files detected. Stored as read1 = ".$$sample{'reads1FqGz'}." read2 = ".$$sample{'reads2FqGz'};
					#warn Dumper($sample);
				#se without index
				}elsif(index($fqbase,"_R1_") && ( -e glob($prefix."/".$fqbase.$ext) && -e glob($prefix."/".$fqbase2.$ext) )){
					
					$sample -> {"reads1FqGz"}=glob($prefix."/".$fqbase.$ext);
					$sample -> {"reads1FqGzMd5"}=Md5Sum($sample -> {"reads1FqGz"});

					$sample -> {"reads2FqGz"}="";
					$sample -> {"reads3FqGz"}=glob($prefix."/".$fqbase2.$ext);
					$sample -> {"reads3FqGzMd5"}=Md5Sum($sample -> {"reads3FqGz"});

					warn "## ".localtime(time())." ## INFO ## Single end read files detected. Stored as read1 = ".$$sample{'reads1FqGz'}." read2 = ".$$sample{'reads2FqGz'};
					#warn Dumper($sample);
				}
				else{
					#die "uncaugth read name ";
				}
					
				if($samplePushed == 0){
					#warn '152'.Dumper($sample);
					if (defined($sample -> {"reads1FqGz"})){
						#warn 'push';
						push(@{$samplesheetNew},$sample);
						$samplePushed++;
					}
				}
			}
		}
		$scount++;
		die "## ".localtime(time())." ## ERR ## Sample is missing or not present in $prefix, plz edit samplesheet/prefix!" if(! $samplePushed);
	}
	#warn Dumper($samplesheetNew,$samplesheet);
	warn "## ".localtime(time())." ## INFO ## Samples detected ".scalar(@{$samplesheetNew})."/".scalar(@{$samplesheet})."\n";
	return $samplesheetNew;
}

sub ConvertSampleSheet {
	my $samplesheet = shift @_;
	my $reformatTemplate = shift @_;
	
	my $newSamplesheet;
	
	for my $sample (@$samplesheet){
		my $newSample;
		for my $key (keys(%$reformatTemplate)){
			if(ref $$reformatTemplate{$key} eq 'ARRAY'){
				my $a = $$reformatTemplate{$key};
				for my $val (@$a){
					$$newSample{ $val } = $$sample{$key};
				}
			}elsif($key eq "#MD5SAMPLE"){
				my $string=Dumper($sample);
				my $cmd = "echo '$string'| md5sum | awk '{print \$1}'";
				my $md5sum = `$cmd`;
				chomp $md5sum;
				#warn "md5sum of sample hash with cmd '$cmd' and md5sum '$md5sum'";
				
				$$newSample{ $$reformatTemplate{$key} } = $md5sum;
			}elsif(defined($$sample{$key})){
				$$newSample{ $$reformatTemplate{$key} } = $$sample{$key};
			}else{
				$$newSample{ $$reformatTemplate{$key} } = "";
			}
		}
		push (@$newSamplesheet, $newSample);
	}
	return $newSamplesheet;
}
sub ConvertSampleSheet_OLD {
	my $samplesheet = shift @_;
	my $reformatTemplate = shift @_;
	
	my $newSamplesheet;
	
	for my $sample (@$samplesheet){
		my $newSample;
		for my $key (keys(%$reformatTemplate)){
			if(ref $$reformatTemplate{$key} eq 'ARRAY'){
				my $a = $$reformatTemplate{$key};
				for my $val (@$a){
					$$newSample{ $val } = $$sample{$key};
				}
			}elsif($key eq "#MD5SAMPLE"){
				my $string=Dumper($sample);
				my $cmd = "echo \"".quotemeta($string)."\"| md5sum | awk '{print \$1}'";
				#warn 'something '.quotemeta($string);
				my $md5sum = `$cmd`;
				chomp $md5sum;
				#warn "md5sum of sample hash with cmd '$cmd' and md5sum '$md5sum'";
				
				$$newSample{ $$reformatTemplate{$key} } = $md5sum;
			}elsif(defined($$sample{$key})){
				$$newSample{ $$reformatTemplate{$key} } = $$sample{$key};
			}else{
				$$newSample{ $$reformatTemplate{$key} } = "";
			}
		}
		$$newSample{ "sequencer" } = "illumina";
		
		push (@$newSamplesheet, $newSample);
	}
	return $newSamplesheet;
}

#seqType,sequencerId,sequencingStartDate,reads3FqGz
sub ReformatHashIlluminaToSampleSheet {
	my %reformatHash = ("line"	=>	"internalId",
		"lane"	=>	"lane",
		"seqType"	=>	["sequencerId","seqType"],
		"sequencer"	=>	"sequencer",
		"Sample_ID"	=>	"externalSampleID",
		"Date"	=>	["sequencingStartDate","Date"],
		"run"	=>	"run",
		"flowcellId"	=>	"flowcellId",
		"Assay"	=>	"barcodeType",
		"Assay"	=>	"Assay",
		"index"	=>	"barcode",
		"index1"	=>	"barcode1",
		"index2"	=>	"barcode2", 
		"Sample_ID"  =>      ["sampleName","controlSampleName"],
		#"Sample_Name"	=>	["sampleName","controlSampleName"],
		"IEMFileVersion"	=>	"IEMFileVersion",
		"Experiment Name"	=>	"project",
		"contact"	=>	"contact",
		"Sample Type"	=>	"applicationType",
		"arrayFile"	=>	"arrayFile",
		"arrayID"	=>	"arrayID",
		"capturingKit"	=>	"capturingKit",
		"prepKit"	=>	["prepKit","samplePrep"],
		"AdapterRead2"	=>	"AdapterRead2",
		"I7_Index_ID"	=>	"I7_Index_ID",
		"Adapter"	=>	"Adapter",
		"I5_Index_ID"	=>	"I5_Index_ID",
		"Chemistry"	=>	"Chemistry",
                "Investigator Name"=>"Investigator_Name",
		"reads3FqGz"=>"reads3FqGz",
		"reads2FqGz"=>"reads2FqGz",
		"reads1FqGz"=>"reads1FqGz",
		"reads3FqGzMd5"=>"reads3FqGzMd5",
		"reads2FqGzMd5"=>"reads2FqGzMd5",
		"reads1FqGzMd5"=>"reads1FqGzMd5",

		);
		#warn ref($reformatHash{"sequencer"});#die Dumper(\%reformatHash);
		return \%reformatHash;
}

sub ReformatHashGccToRnaSeq {
	my %reformatHash = ("#MD5SAMPLE"	=>	"internalId",
		"lane"	=>	"lane",
		"sequencer"	=>	["sequencerId","sequencer","seqType"],
		"Sample"	=>	"externalSampleID",
		"sequencingStartDate"	=>	"sequencingStartDate",
		"run"	=>	"run",
		"flowcell"	=>	"flowcellId",
		"seqType"	=>	"singleOrPairedEnd",
		"barcodeType"	=>	"barcodeType",
		"barcode"	=>	"barcode",
		"barcode1"	=>	"barcode1",
		"barcode2"	=>	"barcode2", 
		"externalSampleID"	=>	["sampleName","controlSampleName"],
		"Merge with"	=>	"Merge with",
		"project"	=>	"project",
		"contact"	=>	"contact",
		"Sample Type"	=>	"applicationType",
		"arrayFile"	=>	"arrayFile",
		"arrayID"	=>	"arrayID",
		"capturingKit"	=>	"capturingKit",
		"prepKit"	=>	["prepKit","samplePrep"],
		"GAF_QC_Name"	=>	"GAF_QC_Name",
		"GAF_QC_Date"	=>	"GAF_QC_Date",
		"GAF_QC_Status"	=>	"GAF_QC_Status",
		"GCC_Analysis"	=>	"GCC_Analysis",
		"Plates in stock - DNA sequencing (whole genome / capturing)"	=>	"Plates in stock - DNA sequencing (whole genome / capturing)",
		"Rejected for processing"	=>	"Rejected for processing",
		"Plates in stock - RNA sequencing"	=>	"Plates in stock - RNA sequencing",
		"Barcode 2"	=>	"Barcode 2",
                "reads3FqGz"=>"reads3FqGz",
		"reads2FqGz"=>"reads2FqGz",
		"reads1FqGz"=>"reads1FqGz",
		"reads3FqGzMd5"=>"reads3FqGzMd5",
		"reads2FqGzMd5"=>"reads2FqGzMd5",
		"reads1FqGzMd5"=>"reads1FqGzMd5",

		);
		#warn ref($reformatHash{"sequencer"});#die Dumper(\%reformatHash);
		return \%reformatHash;
}


sub StatsCollecterSampleSheet{
	
	my $samplesheet = shift @_;
	my $statsNew;
	
	for my $sample (@$samplesheet){
		$statsNew -> { 'records'} ++;
		for my $param (keys(%{$sample})){
			#$statsNew -> { 'parameters'} -> {'counts'} -> {$param}++;
			$statsNew -> { 'parameters'} -> {'counts'} -> {$param} -> {$sample -> {$param}}++;
		}
		#project to fastq more than once
		#sample to fastq more than once
		$statsNew -> { 'fishyparametercombinations'} -> {'fastqtoproject'} -> {$sample -> {'project'}."\t".$sample -> {'reads1FqGz'}."\t".$sample -> {'reads2FqGz'}."\t".$sample -> {'reads3FqGz'}}++;
		$statsNew -> { 'fishyparametercombinations'} -> 
			{'fastqtodiffsamplename'} -> 
			{$sample -> {'reads1FqGz'}."\t".$sample -> {'reads2FqGz'}."\t".$sample -> {'reads3FqGz'}} ->
			{$sample -> {'sampleName'}}++;
	}
	#what should this show

	#Everything strange to the samplesheet XD
	my $fishy = '';
	my $fishycount = 0;
	#Multiple references to the same FQ with different samplename.
	for my $projectFastqCombi (keys(%{$statsNew -> { 'fishyparametercombinations'} -> {'fastqtoproject'}})){
		if($statsNew -> { 'fishyparametercombinations'} -> {'fastqtoproject'} -> { $projectFastqCombi} > 1){
			$fishy .= "$projectFastqCombi\n" if($fishycount < 50);
			$fishycount++;
		}
		
	}
	print "## following project/fastq combinations[max 50 of $fishycount] indicate the use of the same fastq in the same project multiple times please inspect.\nproject\tR1\tR2\tR3\n$fishy" if($fishy);
	
	$fishy = '';
	$fishycount = 0;
	#sampleName to FastqDuplicates
	for my $fastqtodiffsamplename (keys(%{$statsNew -> { 'fishyparametercombinations'} -> {'fastqtodiffsamplename'}})){
		if(scalar(keys(%{$statsNew -> { 'fishyparametercombinations'} -> {'fastqtodiffsamplename'} -> {$fastqtodiffsamplename}})) > 1){
			if(scalar(keys(%{$statsNew -> { 'fishyparametercombinations'} -> {'fastqtodiffsamplename'} -> {$fastqtodiffsamplename}})) > 1 && $fishycount < 50){
				$fishy .=$fastqtodiffsamplename."\t".join("\t",keys(%{$statsNew -> { 'fishyparametercombinations'} -> {'fastqtodiffsamplename'} -> {$fastqtodiffsamplename}}))."\n";
			}
			$fishycount++;
		}
	}
	print "## following fastq/sample combinations[max 50 of $fishycount] indicate inconsistent sample naming possibly between projects.\nR1\tR2\tR3\tsampleName\t[sampleName]\n$fishy" if($fishy);
	
	#warn Dumper($statsNew -> { 'parameters'} -> {'counts'})." ";
	
	warn "## ".localtime(time())." ## INFO ## Samples detected ".scalar(@$samplesheet)."\n";
	return $statsNew;
}
sub ReadSamplesheetWithValidation {
	my $samplesheetf = shift @_;
	my $mods = shift @_;
	my $samplesheet;
	warn "[INFO] Analysing file '$samplesheetf'.\n";
	open( my $samplesheeth,'<', $samplesheetf) or die "[ERROR] Cannot read samplesheet file:'$samplesheetf'";
	$_=<$samplesheeth>;
	chomp;
	my @h = CommaseparatedSplit($_);
	#die Dumper(\@h);
	while(<$samplesheeth>){
		chomp;
		my @c = CommaseparatedSplit($_);
		ValidateColvalues(\@h,\@c) if($mods -> {'cols'});
		ValidateSpecialChars(\@h,\@c) if($mods -> {'chars'});
		#die Dumper(\@c);
		my %d;
		my $i=0;
		map{$d{$_}=$c[$i]; $i++;}(@h);
		$i++;
		#$c[$i]=join(",",@h);
		#ReadFileNameConstructor(\%d);
		push(@$samplesheet,\%d);

	}
	return $samplesheet;
}
sub ValidateColvalues {
	my $header = shift @_;
	ValidateHeader($header);
	my $columns = shift @_;
	if(scalar(@{$header}) ne scalar(@{$columns})){
		die "[VALIDATIONERROR] number of columns (".scalar(@{$columns}).") in line $. are not equal to columns in header ".scalar(@{$columns})
			.".\nArray dump of header".Dumper($header)
			.".\nArray dump of $. columns".Dumper($columns)." ";
	}
}
sub ValidateSpecialChars {
        my $header = shift @_;
        my $columns = shift @_;
	for my $field (@{$header}, @{$columns}){
		if($field =~ /[ \&\%\@\$\'\"\[\]\{\}\*\!\~\`]/){
			die "[VALIDATIONERROR] field '".$field."' matches /[ \#\&\%\@\$\-\'\"\[\]\{\}\*\!\~\`]/ in line '$.' "
                	        .".\nArray dump of header".Dumper($header)
                	        .".\nArray dump of $. columns".Dumper($columns)." ";
		}
                if($field =~ /[\-]/){
                        warn "[VALIDATIONWARNING] field '".$field."' matches /[\-]/ in line '$.'. Consider removing.";
                }

	}
}
sub ValidateHeader {
	my $header = shift(@_);
	
	my @requiredValues = (
		'internalId',
		'samplePrep',
		'seqType',
		'sequencer',
		'sequencerId',
		'run',
		'flowcellId',
		'lane',
		'sampleName',
		'barcode',
		'sequencingStartDate',
		'project',
		'controlSampleName',
		'reads1FqGz',
		'reads3FqGz',
		'reads2FqGz');
	
	my @dieRequired = ();
	for my $requirement (@requiredValues){
		push(@dieRequired,$requirement) if(scalar(grep(${requirement} eq $_,@{$header})) != 1);
		warn "[VALIDATIONWARNING] is this ok? $requirement ".join(',',(@{$header}, scalar(grep(${requirement} eq $_,@{$header})) )) if(scalar(grep(${requirement} eq $_ , @{$header})) != 1);
	}
	die "[VALIDATIONERROR] Missing required field(s) or field(s) declared twice: '".join(',',@dieRequired)."' in samplesheet header '".join(',',@{$header})."'" if(scalar(@dieRequired)>0);
	
	
	#(
	#	'reads1FqGz',
	#	'reads3FqGz',
	#	'reads2FqGz'
	#);
	
	#@{$header}
	#die Dumper($header);
}

sub ValidateSequence{
	#($mergedSamplesheet);#ILLUMINA,SLX,SOLEXA,SOLID,454,LS454,COMPLETE,PACBIO,IONTORRENT,CAPILLARY,HELICOS,UNKNOWN

	my $samplesheet = shift @_;
	
	for my $sample (@{$samplesheet}){
		die "[VALIDATIONERROR] Invalid 'sequencer' field contains '". uc($sample -> {'sequencer'}).
			"' should match ILLUMINA, SLX, SOLEXA, SOLID, 454, LS454, COMPLETE, PACBIO, IONTORRENT, CAPILLARY, HELICOS or UNKNOWN"
			 if(not(uc($sample -> {'sequencer'}) =~ m/ILLUMINA|SLX|SOLEXA|SOLID|454|LS454|COMPLETE|PACBIO|IONTORRENT|CAPILLARY|HELICOS|UNKNOWN/ ));
	}
}

sub ValidateControlSampleNames {
	my $samplesheet = shift @_;
	
	#my $samplesh	
	return if (not(defined($samplesheet -> [0] -> {'controlSampleName'})));
	my %combihash;
	my $s;#sample
	#gather array for checking the controlsample per project and one for checking amount of combinations
	map{push(@{$s}, $_ -> {'project'} . "|" . $_ -> {'sampleName'});
		$combihash{ $_ -> {'project'} }{'samples'}{ $_ -> {'sampleName'} }++;
		$combihash{ $_ -> {'project'} }{'controlsampleNames'}{ $_ -> {'controlSampleName'} }++;
		$combihash{ $_ -> {'project'} }{'samplecontrol'}{ $_ -> {'sampleName'} . '|' . $_ -> {'controlSampleName'} }++;
	}(@{$samplesheet});
	my $i=1;
	my $cs;	#controlsample
	map{push(@{$cs}, {'project' => $_ -> {'project'}, 'controlSampleName' => $_ -> {'controlSampleName'}, 'line' =>  $i}); $i++;}(@{$samplesheet});
	my $error = "";
	map{
		my $csname = $_ -> {'controlSampleName'};
		my $csproject = $_ -> {'project'};
		my $csline = $_ -> {'line'};
		my $seen = 0;
		map{
			$seen++ if($_ eq $csproject . '|' . $csname);
		}(@{$s});
		$error = "'".$csname."' in project '".$csproject."' at $csline," if(not($seen));
	}(@{$cs});
	
	#reduce the amount of normals
	my $maxcontrol = 3;
	for my $project (keys %combihash){
		warn "[VALIDATIONWARNING] project '".$project.
			"' has amount of samples '".scalar(keys(%{$combihash{$project}{'samples'}})).
			"' has amount of controls '".scalar(keys(%{$combihash{$project}{'controlsampleNames'}})).
			"' amount of controlsamples/sample combinations '".scalar(keys(%{$combihash{$project}{'samplecontrol'}}))."'\n";
		die "[VALIDATIONERROR] in project '$project' amount of controlsamples '".scalar(keys(%{$combihash{$project}{'controlsampleName'}}))."' exceeds '$maxcontrol'" if( scalar(keys(%{$combihash{$project}{'controlsampleName'}})) > $maxcontrol); 
	}
	
	chop $error;
	die "[VALIDATIONERROR] The following samplename(s) are seen in the controlsamplenames row but not seen in the samplenames row for the project (format'controlsaplename' at \$lineno) : \n".$error if($error ne "");
	
	
}
sub ValidateFqFiles {
	my $samplesheet = shift @_;
	for my $fileparam ('reads1FqGz',
        	'reads3FqGz',
        	'reads2FqGz'){
		#next if (not(defined($samplesheet -> [0] -> {$file})));
		for my $sample (@{$samplesheet}){
			next if(not(defined($sample -> {$fileparam})));
			die "[VALIDATIONERROR] Invalid file '".$sample -> {$fileparam}."' in parameter '".$fileparam."' ".Dumper($sample) 
				if(defined($sample -> {$fileparam}) &&  $sample -> {$fileparam} ne "" && ! -e $sample -> {$fileparam});
		}
	}
}
sub ValidateInternalId {
	my $samplesheet = shift @_;
	for my $param ('internalId'){
		my %seen;
		#next if (not(defined($samplesheet -> [0] -> {$file})));
		for my $sample (@{$samplesheet}){
			next if(not(defined($sample -> {$param})));
			die "[VALIDATIONERROR] Invalid sampleid/file mapping seen internalId id'".
				$sample -> {$param}.
				"' sample'".$sample -> {'sampleName'}."' more than once with different files: '".
				$sample -> {'reads1FqGz'}.
				"' & '".$seen{ $sample -> {$param} }{ $sample -> {'sampleName'} }{'fileold'}.
				"'" 
				if(defined($seen{ $sample -> {$param} }{ $sample -> {'sampleName'} }{'seen'}) && 
				defined($seen{ $sample -> {$param} }{ $sample -> {'sampleName'} }{'fileold'}) &&
				$seen{ $sample -> {$param} }{ $sample -> {'sampleName'} }{'fileold'} ne $sample -> {'reads1FqGz'});
			$seen{ $sample -> {$param} }{ $sample -> {'sampleName'} }{'seen'}++;
			$seen{ $sample -> {$param} }{ $sample -> {'sampleName'} }{'fileold'}=$sample -> {'reads1FqGz'};
		}
		warn Dumper(\%seen)." ";
	}
}

sub HashWalkBestSamples {
	my ($hash, $key_list, $samplelist, $callback) = @_;
	while (my ($k, $v) = each %$hash) {
		# Keep track of the hierarchy of keys, in case
		# our callback needs it.
		push @$key_list, $k;
	
		if (ref($v) eq 'HASH') {
			# Recurse.
			HashWalkBestSamples($v, $key_list, $samplelist, $callback);
		}
		else {
			# Otherwise, invoke our callback, passing it
			# the current key and value, along with the
			# full parentage of that key.
			$callback->($k, $v, $samplelist, $key_list);
		}
		pop @$key_list;
	}	
}

sub print_keys_and_value {
    my ($k, $v, $key_list) = @_;
    printf "k = %-8s  v = %-4s  key_list = [%s]\n", $k, $v, "@$key_list";
}

sub CollectBestSamples {
	my $key = shift @_;
	my $samples = shift @_;
	my $sampleList = shift @_;
	my $keylist = shift @_;
	my $bestSample= shift @{$samples};
	if(scalar @{$samples}){
		for my $sample (@{$samples}){
			if(length($sample -> { 'sampleName'}) > length($bestSample -> { 'sampleName'}) && index($sample -> { 'sampleName'},$bestSample -> { 'sampleName'}) > 0){
				$bestSample = $sample;
			}elsif($sample -> { 'sampleName'} eq $bestSample -> { 'sampleName'}){
				#keep original if equal
				#die "Dem names are eq!".Dumper($sample,$bestSample) ." ";
			}else{
				#die "Exeptions are there to be caugth".Dumper($sample,$bestSample) . " ";
			}
		}
	}
	push(@{$sampleList}, $bestSample);
}

sub Batcher{
	#could oalso be named botcher cause this can botch up yer samplesheet
	my ($samplesheet, $projectbase, $maxbatchsize, $controlSampleName, $moreControls) = @_;
	warn "Batcher".Dumper(@_). ' ';
	my $controlsample = GetControl($samplesheet, $controlSampleName);
	warn "##### $moreControls";
	my @moreControls = split(',',$moreControls);
	my $batchindex = 0;
	warn "here : ".Dumper($controlsample)."\n";
	my $cs;push(@{$cs},$controlsample);
	
	my $batchedSamplesheet;
	my @batchedProjectNames;
	my %sampleNameToBatch; #for checkin if samplename is already seen and adding samplename to the correct batch
	for my $sample (@{$samplesheet}){
		if(defined($sampleNameToBatch{$sample -> {'sampleName'}})){
			
			$sample -> {'project'} = $sampleNameToBatch{$sample -> {'sampleName'}};
			$sample -> {'controlSampleName'} = $controlsample -> {'sampleName'};
			push(@{$batchedSamplesheet}, $sample);
		}else{
			my $batchProjectName = $projectbase . (floor($batchindex / $maxbatchsize));
			
			$sampleNameToBatch{$sample -> {'sampleName'}} = $batchProjectName;

			if((floor($batchindex / $maxbatchsize) == $batchindex / $maxbatchsize)){
				push(@batchedProjectNames,$batchProjectName);
				warn "INFO: Adding new project $batchProjectName\n";
				my $newControl; %{$newControl} =  %{$controlsample};			
				$newControl -> {'project'} = $batchProjectName;
				$newControl -> {'controlSampleName'} = $newControl -> {'sampleName'};
				push(@{$batchedSamplesheet}, $newControl);
			}
			$sample -> {'project'} = $batchProjectName;
			$sample -> {'controlSampleName'} = $controlsample -> {'sampleName'};
			push(@{$batchedSamplesheet}, $sample);
			#print $projectbase. (floor($batchindex / $maxbatchsize)) . SamplesheetAsString($cs).Dumper($cs)."\n";
			$batchindex++;
		}
	}
	#now cleanup any second 'controlsamples'
	my $sampleindex = 0;
	my $seen;
	while($sampleindex < scalar(@{$batchedSamplesheet})){
		if($batchedSamplesheet -> [$sampleindex ] -> {'sampleName'} eq $controlsample -> {'sampleName'} ){
			$seen -> {$batchedSamplesheet -> [$sampleindex ] -> {'sampleName'}} -> {$batchedSamplesheet -> [$sampleindex ] -> {'project'}}++;
			if($seen -> {$batchedSamplesheet -> [$sampleindex ] -> {'sampleName'}} -> {$batchedSamplesheet -> [$sampleindex ] -> {'project'}} >=2){
				splice(@{$batchedSamplesheet},$sampleindex,1);
			}else{
				$sampleindex++;		
			}
			
		}else{
			$sampleindex++;
		}
	}
	for my $batch (@batchedProjectNames){
		for my $control (@moreControls){
			if($control ne $controlSampleName){
				my $controldat = GetControl($batchedSamplesheet, $control);
				$controldat -> {'controlSampleName'} = $controldat -> {'sampleName'};
				$controldat -> {'project'} = $batch;
				push(@{$batchedSamplesheet}, $controldat);
			}		
		}
	}
	return $batchedSamplesheet;
}
sub GetControl {
	my ($samplesheet,$controlSampleName) = @_;
	for my $sample (@{$samplesheet}){
		#warn "sample ".$sample -> {'sampleName'}.".".$controlSampleName.".";
		#
		if( $sample -> {'sampleName'} eq $controlSampleName){
			#warn "##################################################iseq!!!!!!!!!!!!!!!";
			my $sampleNew; %{$sampleNew}=%{$sample};
			return $sampleNew;
		}
	}
	die "Could not find controlsample '$controlSampleName'.";
	
}
sub SplitByProject {
	my ($samplesheet) = @_;
	my $last;
	my @splitsamplesheet;
	my $byProject;
	for my $sample (sort {$a -> {'project'} cmp $b -> {'project'}} (@{$samplesheet})){
		push(@{$byProject -> {$sample -> {'project'}}},$sample);
		#if(not(defined($last))){
		#	push(@{$splitsamplesheet[0]},$sample)
		#}elsif($last eq $sample -> {'project'}){
		#	push(@{$splitsamplesheet[-1]},$sample)
		#}else{
		#	push(@splitsamplesheet,[$sample]);
		#}
		#
		#$last = $sample -> {'project'};
	}
	return values(%{$byProject});
}

sub CmdRunner {
    my $ret;
    my $cmd = join(" ",@_);

    warn localtime( time() ). " [INFO] system call:'". $cmd."'.\n";

    @{$ret} = `($cmd )2>&1`;
    if ($? == -1) {
        die localtime( time() ). " [ERROR] failed to execute: $!\n";
    }elsif ($? & 127) {
        die localtime( time() ). " [ERROR] " .sprintf "child died with signal %d, %s coredump",
         ($? & 127),  ($? & 128) ? 'with' : 'without';
    }elsif ($? != 0) {
        die localtime( time() ). " [ERROR] " .sprintf "child died with signal %d, %s coredump",
             ($? & 127),  ($? & 128) ? 'with' : 'without';
    }else {
        warn localtime( time() ). " [INFO] " . sprintf "child exited with value %d\n", $? >> 8;
    }
    return @{$ret};
}

sub Md5Sum {
	my $file = $_[0];
	my $ret;@{$ret} = CmdRunner("md5sum '$file'");
	return substr(${$ret}[0],0,32);
}
