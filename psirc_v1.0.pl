#!/usr/bin/perl
use strict; use warnings;
use Getopt::Long qw(GetOptions);

##
#BEGIN { unshift @INC, "/home/hoyu/scripts"; }
#use DataBrowser;
##

my $version = "_v1.0";
my $kallisto_fork = "kallisto_b";

#Do these sam flags indicate they are first or second in pair?
my %in_pair = (
	99 => '1',
	355 => '1',
	147 => '2',
	403 => '2',
	83 => '1',
	339 => '1',
	163 => '2',
	419 => '2',
	69 => '1',
	73 => '1',
	89 => '1',
	101 => '1',
	329 => '1',
	345 => '1',
	77 => '1',
	133 => '2',
	137 => '2',
	153 => '2',
	165 => '2',
	393 => '2',
	409 => '2',
	141 => '2',
);

my %transcripts_seen;
my %bsjs_with_equiv_fsj;

my $kallisto_k = 31;
my $kallisto_l = 30;
my $segment_kallisto_l = 15;

my $gen_indices;
my $threads_kallisto = 1;
my $EXPECT_READIDS_PER_BSJ = 2;
#my $DEFAULT_UNIQUE_BY = 10;
my $MBSL = 7;
my $output_dir = ".";
my $final_kallisto_k = 21;
my $straight_to_full_length;
my $equiv_fsj;
my $lib_type = "fr-unstranded";
my $unique_by = 10;
my $linear_ratio_info;
my $max_mappable_bsjs = 300;
my $max_mappable_alt_fsjs = 20;
my $max_alt_fsj_per_isoform = 5;
my $full_length;

GetOptions(
	'i|indices' => \$gen_indices,
	't|threads=i' => \$threads_kallisto,
	'o|output-dir=s' => \$output_dir,
	'y|library-type=s' => \$lib_type,
	'm|mbsl=i' => \$MBSL,
	'r|final-k=i' => \$final_kallisto_k,
	'u|unique-by=i' => \$unique_by,
	's|straight-to-full-length' => \$straight_to_full_length,
	'j|equiv-fsj' => \$equiv_fsj,
	'l|linear-ratio-info' => \$linear_ratio_info,
	'b|max-mappable-bsjs=i' => \$max_mappable_bsjs,
	'a|max-mappable-alt-fsjs=i' => \$max_mappable_alt_fsjs,
	'p|max-alt-fsj-per-isoform=i' => \$max_alt_fsj_per_isoform,
	'f|full-length' => \$full_length,
) or die_usage();

#segment partition settings
my ($p_1_kmer, $p_1_kmer_bp_low, $p_1_kmer_bp_high, $p_1_kmer_target_length) = ($MBSL, $MBSL, 14, 20);
my ($p_2_kmer, $p_2_kmer_bp_low, $p_2_kmer_bp_high, $p_2_kmer_target_length) = (13, 15, 22, 30);
my ($p_3_kmer, $p_3_kmer_bp_low) = (21, 23);

die "$output_dir does not exist\n" unless -d "$output_dir";

die "Supported library types are fr-unstranded, fr-firststrand, or fr-secondstrand\n" 
unless ( ($lib_type eq "fr-unstranded") || ($lib_type eq "fr-firststrand") || ($lib_type eq "fr-secondstrand") );

die "-m/--mbsl must be a whole number between 5 and 10\n" unless ( ($MBSL >= 5) && ($MBSL <= 10) );

die "-r/--final-k must be an odd whole number between 1 and 31\n" unless ( ($final_kallisto_k % 2 == 1) 
&& ($final_kallisto_k >= 1) && ($final_kallisto_k <= 31) && ($final_kallisto_k =~ /^\d+\z/) );

die "-f/--full-length cannot be used with -a/--linear-ratio-info. Use the standard run or -s instead\n"
if ($linear_ratio_info && $full_length);

my ($fa, $kallisto, $fq1, $fq2, $provided_dir);

if ($gen_indices) {
	die_usage_4() unless (scalar(@ARGV)==2);
	die "fasta file names need to end with .fa or .fasta\n" unless ($ARGV[0] =~ /.*(\.fa|\.fasta)$/);
	die "$ARGV[0] not found\n" unless (-e $ARGV[0]);
	die "$ARGV[1] not found\n" unless (-e $ARGV[1]);
	($fa, $kallisto) = @ARGV;
} elsif ($linear_ratio_info) {
	die_usage_3() unless (scalar(@ARGV)==5);
	die "fasta file names need to end with .fa or .fasta\n" unless ($ARGV[1] =~ /.*(\.fa|\.fasta)$/);
	
	die "$ARGV[0] not found\n" unless (-e $ARGV[0]);
	die "$ARGV[1] not found\n" unless (-e $ARGV[1]);
	die "$ARGV[2] not found\n" unless (-e $ARGV[2]);
	die "$ARGV[3] not found\n" unless (-e $ARGV[3]);
	die "$ARGV[4] not found\n" unless (-e $ARGV[4]);
	($provided_dir, $fa, $kallisto, $fq1, $fq2) = @ARGV;
} elsif ($straight_to_full_length) {
	die_usage_2() unless (scalar(@ARGV)==5);
	die "fasta file names need to end with .fa or .fasta\n" unless ($ARGV[1] =~ /.*(\.fa|\.fasta)$/);
	
	die "$ARGV[0] not found\n" unless (-e $ARGV[0]);
	die "$ARGV[1] not found\n" unless (-e $ARGV[1]);
	die "$ARGV[2] not found\n" unless (-e $ARGV[2]);
	die "$ARGV[3] not found\n" unless (-e $ARGV[3]);
	die "$ARGV[4] not found\n" unless (-e $ARGV[4]);
	
	($provided_dir, $fa, $kallisto, $fq1, $fq2) = @ARGV;
}
 else {
	die_usage() unless (scalar(@ARGV)==4);
	
	die "$ARGV[0] not found\n" unless (-e $ARGV[0]);
	die "fasta file names need to end with .fa or .fasta\n" unless ($ARGV[0] =~ /.*(\.fa|\.fasta)$/);
	die "$ARGV[1] not found\n" unless (-e $ARGV[1]);
	die "$ARGV[2] not found\n" unless (-e $ARGV[2]);
	die "$ARGV[3] not found\n" unless (-e $ARGV[3]);

	($fa, $kallisto, $fq1, $fq2) = @ARGV;
}

if ($straight_to_full_length || $linear_ratio_info) {
	$output_dir = $provided_dir;
}

#Remove final slash in output_dir if given
if ($output_dir =~ m/\/$/) {
	chop($output_dir);
}

my $mapping_dir = $output_dir . "/mapping";
unless ($gen_indices) {
	mkdir $mapping_dir unless (-d $mapping_dir);
}

my $target_potential_backsplice_fasta = $mapping_dir . "/" . "potential_backsplice_transcripts.fa";
my $target_potential_backsplice_fasta_no_copied_seq = $mapping_dir . "/" . "potential_backsplice_transcripts.no_copied_seq.fa";
my $potential_backsplice_reads_1 = $mapping_dir . "/possible_bsj_reads_1.fastq.gz";
my $potential_backsplice_reads_2 = $mapping_dir . "/possible_bsj_reads_2.fastq.gz";

my $low_num_temp_sam = $mapping_dir . "/temp_final_mapping.sam";
my $num_transcripts;
my $low_num_transcripts = 500;
my $num_reads;
my $low_num_reads = 100000;
my $low_num;
my $low_num_sam=0;
my %low_num_read_ids;

my ($fa_file_pre, $full_path_fa_file_pre);
if ($fa =~ m/(.*\/)(.*\.fasta)$/) {
	my $path = $1;
	my $fa_file = $2;
	$fa_file =~ s/\.fasta//;
	$fa_file_pre = $fa_file;
	$full_path_fa_file_pre = $path . $fa_file_pre;
} elsif ($fa =~ m/(.*\/)(.*\.fa)$/) {
	my $path = $1;
	my $fa_file = $2;
	$fa_file =~ s/\.fa//;
	$fa_file_pre = $fa_file;
	$full_path_fa_file_pre = $path . $fa_file_pre;
} elsif ($fa =~ m/(.*\.fasta)$/) {
	my $fa_file = $1;
	$fa_file =~ s/\.fasta//;
	$fa_file_pre = $fa_file;
	$full_path_fa_file_pre = $fa_file_pre;
} elsif ($fa =~ m/(.*\.fa)$/) {
	my $fa_file = $1;
	$fa_file =~ s/\.fa//;
	$fa_file_pre = $fa_file;
	$full_path_fa_file_pre = $fa_file_pre;
} 

my $all_possible_bsj_targets_fa = "$full_path_fa_file_pre.all_possible_bsj_targets.fa";
my $FirstandLastExons_entities_fa = "$full_path_fa_file_pre.FirstandLastExons_entities.fa";

my $all_possible_bsj_targets_fa_index = $all_possible_bsj_targets_fa . ".index";
my $FirstandLastExons_entities_fa_index = $FirstandLastExons_entities_fa . ".index";

my $READ_LENGTHS;

if ( ($gen_indices) && (!$linear_ratio_info) && (!$straight_to_full_length)) {
	print "Generating all possible bsj and FirstandLastExons entities fa files and their kallisto indices... ";
	system("date");
	generate_targets_and_indices();
	system("date");
	die "Finished generating. The indices should be kept in the same file path as the <transcriptome_fa> and their names should not be changed.\n";
} elsif ( ($gen_indices) && ($linear_ratio_info) ) {
	die_usage_4();
}

if ( (!$linear_ratio_info) && (!$straight_to_full_length) ) {
	$READ_LENGTHS = map_reads_to_targets();
} elsif ( ($straight_to_full_length) && (!$linear_ratio_info) ) {
	print "Getting the required variables... ";
	system("date");
	get_required_variables("s");
	print "Ready. ";
	system("date");
} elsif ( (!$straight_to_full_length) && ($linear_ratio_info) ) {
	print "Getting the required variables... ";
	system("date");
	get_required_variables("l");
	print "Ready. ";
	system("date");
}


print "Final mapping... ";
system("date");


###
if ($low_num) {
	if ($low_num_sam != 1) {
		my $backsplice_transcripts_fa_index = $target_potential_backsplice_fasta . ".index";
		system("$kallisto index -k $final_kallisto_k -i $backsplice_transcripts_fa_index $target_potential_backsplice_fasta");
		system("$kallisto pseudo -i $backsplice_transcripts_fa_index -o $mapping_dir -t $threads_kallisto --pseudobam $potential_backsplice_reads_1 $potential_backsplice_reads_2 > $low_num_temp_sam");
		check_temp_sam($low_num_temp_sam);
	}
}
###


my %transcripts_bsjs;
my %outputted_raw_names;
my %raw_name_strand;
final_mapping_and_output_candidate_circRNAs() unless ($straight_to_full_length);
print "Finished mapping. ";
system("date");
print "\n";



#full length
my $full_length_dir = $output_dir . "/full_length";
my $fsj_targets_fa = $full_length_dir . "/$fa_file_pre" . ".within_backspliced.fsj_targets.fa";
my $target_potential_fli_fasta = $full_length_dir . "/" . "potential_full_length_isoforms.fa";
my $invalidated_flis = $full_length_dir . "/". "invalidated_full_length_isoforms.txt";
my %fli_breakpoints;
my %fli_strands;
my ($potential_fsj_reads_1, $potential_fsj_reads_2) = ("$full_length_dir/potential_fsj_reads_1.fastq.gz", "$full_length_dir/potential_fsj_reads_2.fastq.gz");

$low_num_temp_sam = $full_length_dir . "/temp_final_mapping.sam";
$low_num="";
$low_num_sam=0;


if ($full_length || $straight_to_full_length) {
	print "Full length isoform inference... ";
	system("date");
	mkdir $full_length_dir unless (-d $full_length_dir);

	generate_alt_fsj_targets();
	map_reads_to_targets_fsj();
	
	if ($low_num) {
		my $full_length_isoforms_fa_index = $target_potential_fli_fasta . ".index";
		system("$kallisto index -k $final_kallisto_k -i $full_length_isoforms_fa_index $target_potential_fli_fasta");
		system("$kallisto pseudo -i $full_length_isoforms_fa_index -o $full_length_dir -t $threads_kallisto --pseudobam $potential_fsj_reads_1 $potential_fsj_reads_2 > $low_num_temp_sam");
		check_temp_sam($low_num_temp_sam);
	}
	
	print "Final mapping... ";
	system("date");
	final_mapping_and_output_fsjs();
	print "Finished full length isoform inference. ";
	system("date");
} else {
	exit;
}

sub die_usage_4 {
die "
Usage: perl psirc$version.pl -i/--indices <transcriptome_fa> <$kallisto_fork" . "_executable>    generates the two kallisto indices  
                                                                                     required for psirc input
                                                                                     [outputs to same location as <transcriptome_fa>]
";
}

sub die_usage {
die "
Usage: perl psirc$version.pl [Options] <transcriptome_fa> <$kallisto_fork" . "_executable> <first_in_pair_fastq> <second_in_pair_fastq>

Before using the (above) main script, run the following once to generate the required kallisto indices for <transcriptome_fa>,
Generate indices: perl psirc$version.pl -i/--indices <transcriptome_fa> <$kallisto_fork" . "_executable>    generates the two kallisto indices
                                                                                                required for psirc input
                                                                                                [outputs to same location as <transcriptome_fa>]                        

Options:
-o/--output-dir	<string>        the output directory                                                      [default ./]
-t/--threads <int>              number of threads to run kallisto with                                    [default $threads_kallisto]
-u/--unique-by <int>            max number of genomic bsj loci a read is allowed to be mapped to          [default $unique_by]
-y/--library-type <lib_type>    the library type of the paired input reads                                [default $lib_type]
-m/--mbsl <int>                 minimum back-splice segment length                                        [default $MBSL]
-r/--final-k <int>              the kallisto k-mer size used for final remapping                          [default $final_kallisto_k]
-j/--equiv-fsj                  check if the output the candidate circRNA junction also supports an fsj   [default FALSE]
-s/--straight-to-full-length    straight to full-length isoform inference by providing default run outdir [default FALSE]
-l/--linear-ratio-info          use all the original reads in final step, output linear read counts file  [default FALSE]
-b/--max-mappable-bsjs          maximum number of mappable back-splice junctions for a transcript         [default $max_mappable_bsjs]
-a/--max-mappable-alt-fsjs      maximum number of mappable alt forward-splice junctions for a transcript  [default $max_mappable_alt_fsjs]
-p/--max-alt-fsj-per-isoform    maximum alt forward-splice junctions allowed per transcript isoform       [default $max_alt_fsj_per_isoform]
-f/--full-length                infer full-length circRNA isoforms **necessary for input to psirc-quant** [default FALSE]

<lib_type>:
fr-unstranded
fr-firststrand
fr-secondstrand
";
}

sub die_usage_2 {
die "
Usage: perl psirc$version.pl -s/--straight-to-full-length [Options] <default_run_output_dir> <transcriptome_fa> <$kallisto_fork" . "_executable> <first_in_pair_fastq> <second_in_pair_fastq>

Options:
-t/--threads <int>              number of threads to run kallisto with                                    [default $threads_kallisto]
-u/--unique-by <int>            max number of genomic bsj loci a read is allowed to be mapped to          [default $unique_by]
-y/--library-type <lib_type>    the library type of the paired input reads                                [default $lib_type]
-m/--mbsl <int>                 minimum back-splice segment length                                        [default $MBSL]
-r/--final-k <int>              the kallisto k-mer size used for final remapping                          [default $final_kallisto_k]
-j/--equiv-fsj                  check if the output the candidate circRNA junction also supports an fsj   [default FALSE]
-l/--linear-ratio-info          use all the original reads in final step, output linear read counts file  [default FALSE]
-b/--max-mappable-bsjs          maximum number of mappable back-splice junctions for a transcript         [default $max_mappable_bsjs]
-a/--max-mappable-alt-fsjs      maximum number of mappable alt forward-splice junctions for a transcript  [default $max_mappable_alt_fsjs]
-p/--max-alt-fsj-per-isoform    maximum alt forward-splice junctions allowed per transcript isoform       [default $max_alt_fsj_per_isoform]

<lib_type>:
fr-unstranded
fr-firststrand
fr-secondstrand
";
}

sub die_usage_3 {
die "
Usage: perl psirc$version.pl -l/--linear-ratio-info [Options] <default_run_output_dir> <transcriptome_fa> <$kallisto_fork" . "_executable> <first_in_pair_fastq> <second_in_pair_fastq>

Options:
-t/--threads <int>              number of threads to run kallisto with                                    [default $threads_kallisto]
-u/--unique-by <int>            max number of genomic bsj loci a read is allowed to be mapped to          [default $unique_by]
-y/--library-type <lib_type>    the library type of the paired input reads                                [default $lib_type]
-m/--mbsl <int>                 minimum back-splice segment length                                        [default $MBSL]
-r/--final-k <int>              the kallisto k-mer size used for final remapping                          [default $final_kallisto_k]
-j/--equiv-fsj                  check if the output the candidate circRNA junction also supports an fsj   [default FALSE]
-b/--max-mappable-bsjs          maximum number of mappable back-splice junctions for a transcript         [default $max_mappable_bsjs]
-a/--max-mappable-alt-fsjs      maximum number of mappable alt forward-splice junctions for a transcript  [default $max_mappable_alt_fsjs]
-p/--max-alt-fsj-per-isoform    maximum alt forward-splice junctions allowed per transcript isoform       [default $max_alt_fsj_per_isoform]

<lib_type>:
fr-unstranded
fr-firststrand
fr-secondstrand
";
}


sub get_required_variables {
	my $type=shift;
	#get the hash for existing transcripts first
	open IN, "<$fa" or die "can't open $fa\n";
	while (<IN>) {
		chomp;
		if ($_ =~ m/^.*:\d+-\d+_(.*)_Exon_Lengths_.*_Offsets_.*_/) {
			my $pure_name = $1;
			my ($transcript_name) = $_ =~ m/^>(.*)/;

			if ( (exists $transcripts_seen{$pure_name}) && ($transcripts_seen{$pure_name} ne $transcript_name) ) {
				$transcripts_seen{$transcript_name} = $transcript_name;
				$pure_name = $transcript_name;
			}
			else {
				$transcripts_seen{$pure_name} = $transcript_name;
			}
		}
		
	}
	close IN;
	
	my $outputted_raw_names_list = $mapping_dir . "/outputted_raw_names.lst";
	if ($type eq "s") {
		open IN, "<$outputted_raw_names_list" or die "can't open $outputted_raw_names_list\n";
		while (<IN>) {
			chomp;
			my @line = split("\t", $_);
			my ($raw_name, $num_reads) = ($line[0], $line[1]);
			$outputted_raw_names{$raw_name} = $num_reads;
			my ($pure_name, $bsj);
			if ($raw_name =~ m/^(.*)_(E\d+B\d+)$/) {
				$pure_name = $1;
				$bsj = $2;
			}
			my ($e, $b);
			if ($bsj =~ m/E(\d+)B(\d+)/) {
				$e = $1;
				$b = $2;
			}
			my $transcript_name = $transcripts_seen{$pure_name};
			my $strand;
			if ($transcript_name =~ m/^.*:\d+-\d+_.*_Exon_Lengths_.*_Offsets_.*_(\+|-)/) {
				$strand = $1;
			}
	
			$raw_name_strand{$raw_name} = $strand;

			if ( ($e - $b) == 0 ) {
				next;
			}

			$transcripts_bsjs{$transcript_name}{$bsj}++;	
		}
		close IN;
	}
	my %read_lengths;
	incr_read_lengths(\%read_lengths, "$output_dir/mapping/possible_bsj_reads_1.fastq.gz");
	incr_read_lengths(\%read_lengths, "$output_dir/mapping/possible_bsj_reads_2.fastq.gz");
	
	for my $read_length (sort { $read_lengths{$b} <=> $read_lengths{$a} } keys %read_lengths) {
		$READ_LENGTHS = $read_length;
		last;
	}
	
	print "Most common read length: $READ_LENGTHS\n";
	
}

sub incr_read_lengths {
	my ($read_lengths, $bsj_reads_fq) = @_;
	
	open IN, "gunzip -c $bsj_reads_fq |" or die "can't open $bsj_reads_fq\n";
	my $line=1;
	while (<IN>) {
		if ($line == 2) {
			chomp($_);
			$read_lengths->{length($_)}++;
		}
		if ($line==4) {
			$line=0;
		}
		$line++;
	}
	close IN;
}

sub populate_and_extract_final_step_exon_ends_fasta {
	my ($backsplice_candidate_name_breakpoint_bp, $type) = @_;
	
	my %exon_ends_fasta_p1;
	my %exon_ends_fasta_p2;
	my %exon_ends_fasta_p3;
	
	if ($type eq "bsj") {
		open IN, "<$target_potential_backsplice_fasta" or die "can't open $target_potential_backsplice_fasta\n";
	} elsif ($type eq "fsj") {
		open IN, "<$target_potential_fli_fasta" or die "can't open $target_potential_fli_fasta\n";
	}
	while (<IN>) {
		chomp;
		my $backsplice_candidate_name;
		if ($_ =~ /^>(.*)/) {
			next unless (exists $backsplice_candidate_name_breakpoint_bp->{$1});
			$backsplice_candidate_name = $1;
			my $backsplice_candidate_seq = <IN>;
			
			my @breakpoint_bps;
			my @name_and_begs;
			my @name_and_ends;
			if ($type eq "bsj") {
				push @breakpoint_bps, $backsplice_candidate_name_breakpoint_bp->{$backsplice_candidate_name};
				push @name_and_begs, $backsplice_candidate_name . "&beg";
				push @name_and_ends, $backsplice_candidate_name . "&end";
			} else {
				for my $fsj (keys %{$backsplice_candidate_name_breakpoint_bp->{$backsplice_candidate_name}}) {
					push @breakpoint_bps, $backsplice_candidate_name_breakpoint_bp->{$backsplice_candidate_name}->{$fsj};
					my $name_and_beg = $backsplice_candidate_name . "_" . $fsj . "&beg";
					my $name_and_end = $backsplice_candidate_name . "_" . $fsj . "&end";;
					push @name_and_begs, $name_and_beg;
					push @name_and_ends, $name_and_end;
				}
			}
			
			for (my $i=0; $i < @breakpoint_bps; $i++) {
			
				my ($max_end_seq_offset, $max_end_seq_length, $max_beg_seq_length);
				if ($type eq "bsj") {
					$max_end_seq_offset = 0;
					$max_end_seq_length = $breakpoint_bps[$i];
					$max_beg_seq_length = $breakpoint_bps[$i];
				} elsif ( ($breakpoint_bps[$i] - ($READ_LENGTHS-1)) < 0 ) {
					$max_end_seq_offset = 0;
					$max_end_seq_length = $breakpoint_bps[$i];
					$max_beg_seq_length = $breakpoint_bps[$i];
				} else {
					$max_end_seq_offset = $breakpoint_bps[$i] - ($READ_LENGTHS-1);
					$max_end_seq_length = $READ_LENGTHS-1;
					$max_beg_seq_length = $READ_LENGTHS-1;
				}
			
				my $max_end_seq = substr($backsplice_candidate_seq, $max_end_seq_offset, $max_end_seq_length);
				my $max_beg_seq = substr($backsplice_candidate_seq, $breakpoint_bps[$i], $max_beg_seq_length);
			
				my ($p_1_end_seq, $p_1_beg_seq);
				if ($max_beg_seq_length < $p_1_kmer_target_length) {
					$p_1_end_seq = $max_end_seq;
					$p_1_beg_seq = $max_beg_seq;
				} else {
					$p_1_end_seq = substr($max_end_seq, -($p_1_kmer_target_length));
					$p_1_beg_seq = substr($max_beg_seq, 0, $p_1_kmer_target_length);
				}
			
				my ($p_2_end_seq, $p_2_beg_seq);
				if ($max_beg_seq_length < $p_2_kmer_target_length) {
					$p_2_end_seq = $max_end_seq;
					$p_2_beg_seq = $max_beg_seq
				} else {
					$p_2_end_seq = substr($max_end_seq, -($p_2_kmer_target_length));
					$p_2_beg_seq = substr($max_beg_seq, 0, $p_2_kmer_target_length);
				}
			
				$exon_ends_fasta_p3{$max_end_seq}{$name_and_ends[$i]}++;
				$exon_ends_fasta_p3{$max_beg_seq}{$name_and_begs[$i]}++;
			
				$exon_ends_fasta_p2{$p_2_end_seq}{$name_and_ends[$i]}++;
				$exon_ends_fasta_p2{$p_2_beg_seq}{$name_and_begs[$i]}++;
			
				$exon_ends_fasta_p1{$p_1_end_seq}{$name_and_ends[$i]}++;
				$exon_ends_fasta_p1{$p_1_beg_seq}{$name_and_begs[$i]}++;
			}
		}
	}
	close IN;
	my $output_exon_ends_p1 = extract_final_step_exon_ends_fasta(\%exon_ends_fasta_p1, "1", $type);
	my $output_exon_ends_p2 = extract_final_step_exon_ends_fasta(\%exon_ends_fasta_p2, "2", $type);
	my $output_exon_ends_p3 = extract_final_step_exon_ends_fasta(\%exon_ends_fasta_p3, "3", $type);
	
	return ($output_exon_ends_p1, $output_exon_ends_p2, $output_exon_ends_p3);
}

sub extract_final_step_exon_ends_fasta {
	my ($exon_ends_fasta, $p, $type) = @_;
	my $output_file;
	if ($type eq "bsj") {
		$output_file = $mapping_dir . "/potential_backsplice_transcripts.final_step.junctions.p$p.fa";
	} elsif ($type eq "fsj") {
		$output_file = $full_length_dir . "/fsj.final_step.junctions.p$p.fa";
	}
	open OUT, ">$output_file" or die "can't open $output_file\n";
	for my $beg_or_end_seq (keys %{$exon_ends_fasta}) {
		my @equiv_headers = keys (%{$exon_ends_fasta->{$beg_or_end_seq}});
		my $header = ">" . join("=", @equiv_headers);
		print OUT "$header\n$beg_or_end_seq\n";
	}
	close OUT;
	
	return $output_file;
}

sub print_segment_fastq_seq {

	my ($segment_fastq_seq, $type) = @_;
	
	my $final_remapping;
	if ($type eq "bsj") {
		$final_remapping = $mapping_dir;
	} elsif ($type eq "fsj") {
		$final_remapping = $full_length_dir;
	}

	my $prefix = "$final_remapping/segments_final_step_junctions_";

	my $output_p1_seg = $prefix . $p_1_kmer_bp_low . "to" . $p_1_kmer_bp_high . ".fastq.gz";
	my $output_p2_seg = $prefix . $p_2_kmer_bp_low . "to" . $p_2_kmer_bp_high . ".fastq.gz";
	my $output_p3_seg = $prefix . $p_3_kmer_bp_low . "orabove.fastq.gz";

	open P1, "| gzip -c > $output_p1_seg" or die "can't open $output_p1_seg\n";
	open P2, "| gzip -c > $output_p2_seg" or die "can't open $output_p2_seg\n";
	open P3, "| gzip -c > $output_p3_seg" or die "can't open $output_p3_seg\n";
	
	
	for my $line_2_4 (keys %{$segment_fastq_seq}) {
		for my $read_ids (keys %{$segment_fastq_seq->{$line_2_4}}) {
			my @final_segment_info_ar;
			for my $segment_info (keys %{$segment_fastq_seq->{$line_2_4}->{$read_ids}}) {
				push @final_segment_info_ar, $segment_info;
			}
			my $final_segment_info = join("=", @final_segment_info_ar);
			my $line_1 = $read_ids . "=>" . $final_segment_info;
			
			my @line_2_4_info = split("\n", $line_2_4);
			my $seg_length = length($line_2_4_info[0]);
			if ($seg_length < $p_2_kmer_bp_low) {
				print P1 "@" . $line_1 . "\n" . $line_2_4 . "\n";
			} elsif ($seg_length < $p_3_kmer_bp_low) {
				print P2 "@" . $line_1 . "\n" . $line_2_4 . "\n";
			} else {
				print P3 "@" . $line_1 . "\n" . $line_2_4 . "\n";
			}
			

		}
	}
	close P1;
	close P2;
	close P3;
	
	return ($output_p1_seg, $output_p2_seg, $output_p3_seg);
	
}

sub map_and_examine_final_segments {
	my ($output_exon_ends_fasta, $segments_fastq_file, $segment_kallisto_k, $checked_junctions, $type) = @_;
	
	my $final_remapping;
	if ($type eq "bsj") {
		$final_remapping = $mapping_dir;
	} elsif ($type eq "fsj") {
		$final_remapping = $full_length_dir;
	}
	
	my $output_exon_ends_fa_index = $output_exon_ends_fasta . ".index";
	system("$kallisto index -i $output_exon_ends_fa_index -k $segment_kallisto_k $output_exon_ends_fasta");
	open IN, "$kallisto pseudo -i $output_exon_ends_fa_index -o $final_remapping -t $threads_kallisto --single -l $segment_kallisto_l -s 20 --pseudobam $segments_fastq_file |" or die "can't run kallisto pseudo segments\n";
	my $first_line = <IN>;
	my %target_lengths;
	
	while (<IN>) {
		next if ($_ =~ /^\@PG/);

		if ($_ =~ /^@/) {
			my @line = split("\t", $_);
			chomp(@line);
			my ($length, $target);
			if ($line[2] =~ /LN:(\d+)/) {
				$length = $1;
			}
			 if ($line[1] =~ /SN:(.*)/) {
				$target = $1;
			}
			$target_lengths{$target} = $length;
			next;
		}
		my @line = split(/\s+/, $_);
		my $sam_flag = $line[1];
		next if ($sam_flag == 4);
		my ($segments, $exon_beg_or_ends, $sam_offset, $sam_cigar) = ($line[0], $line[2], $line[3], $line[5]);
		my ($cigar_M) = $sam_cigar =~ /(\d+)M/;
		my @segments_info_ar = split("=>", $segments);
		my $read_id = $segments_info_ar[0];
		my $mapped_segments = $segments_info_ar[1];
		my @mapped_segments_ar = split("=", $mapped_segments);
		
		my @exon_beg_or_ends_ar = split("=", $exon_beg_or_ends);
		
		my @read_id_info = split("&", $read_id);
		my ($in_pair, $raw_read_id) = ($read_id_info[0], $read_id_info[1]);
		
		for my $mapped_segment (@mapped_segments_ar) {
			for my $exon_beg_or_end (@exon_beg_or_ends_ar) {
				if ($mapped_segment eq  $exon_beg_or_end) {
					my @indiv_segment_info_ar = split("&", $mapped_segment);
					my ($backsplice_candidate_name, $segment_beg_or_end) = ($indiv_segment_info_ar[0], $indiv_segment_info_ar[1]);
					##filter out false positives
					if ( ($segment_beg_or_end eq "beg") && ($sam_offset != 1) ) {
						next;
					}
					
					my $target_length = $target_lengths{$exon_beg_or_ends};
					if ($segment_beg_or_end eq "end") {
						if ( ($sam_offset + $cigar_M) != ($target_length+1) ) {
							next;
						}
					}
					
					my $in_pair_beg_or_end = $in_pair . $segment_beg_or_end;
					
					if ($type eq "bsj") {
						$checked_junctions->{$raw_read_id}->{$backsplice_candidate_name}->{$in_pair_beg_or_end}++;
					} else {
						$checked_junctions->{$backsplice_candidate_name}->{$raw_read_id}->{$in_pair_beg_or_end}++;
					}
				}
			}
		}
	}
	close IN;
}

sub final_mapping_and_output_candidate_circRNAs {
	my $backsplice_transcripts_fa_index = $target_potential_backsplice_fasta . ".index";

	system("$kallisto index -k $final_kallisto_k -i $backsplice_transcripts_fa_index $target_potential_backsplice_fasta") unless ($linear_ratio_info);

	my $final_remapping;
	if ($linear_ratio_info) {
		mkdir "$output_dir/linear_ratio_out" unless (-d "$output_dir/linear_ratio_out");
		mkdir "$output_dir/linear_ratio_out/mapping" unless (-d "$output_dir/linear_ratio_out/mapping");
		$output_dir = "$output_dir/linear_ratio_out";
		$final_remapping = "$output_dir/mapping";
	} else {
		$final_remapping = $mapping_dir;
	}
	
	my $proper_pair_out = $final_remapping . "/final_reads_mapped_proper_pair.sam";

	open OUT_P, ">$proper_pair_out" or die "can't open $proper_pair_out\n";
	
	######
	my $sam_headers = $final_remapping . "/sam_headers";
	open OUT_SH, ">$sam_headers" or die "can't open $sam_headers\n";
	######
	
	
	if($linear_ratio_info) {
		open IN, "$kallisto pseudo -i $backsplice_transcripts_fa_index -o $final_remapping -t $threads_kallisto --pseudobam $fq1 $fq2 |" or die "can't run kallisto pseudo\n";
	} else {
		if ( ($low_num) && ($low_num_sam == 1)) {
			open IN, "<$low_num_temp_sam" or die "can't open $low_num_temp_sam\n";
		} elsif ($low_num) {
			open IN, "$kallisto pseudo -i $backsplice_transcripts_fa_index -o $final_remapping -t 1 --pseudobam $potential_backsplice_reads_1 $potential_backsplice_reads_2 |" or die "can't run kallisto pseudo\n";
		} else {
			open IN, "$kallisto pseudo -i $backsplice_transcripts_fa_index -o $final_remapping -t $threads_kallisto --pseudobam $potential_backsplice_reads_1 $potential_backsplice_reads_2 |" or die "can't run kallisto pseudo\n";
		}
		
	}

	
	my %segment_fastq_seq;
			
	my %loci_reads;
	my %loci_tx_exons_names;

	
	my %name_reads;
	my %raw_to_output;
	my %raw_to_loci;
	
	my %potential_olf_reads;
	
	#used for extracting final step exon_ends_fasta
	my %backsplice_candidate_name_breakpoint_bp;
	
	my ($first_strand_expectation, $second_strand_expectation);
	if ($lib_type eq "fr-firststrand") {
		$first_strand_expectation = "-";
		$second_strand_expectation = "+";
	} elsif ($lib_type eq "fr-secondstrand") {
		$first_strand_expectation = "+";
		$second_strand_expectation = "-";
	}
	
	my %linear_name_reads;
	my %linear_loci_reads;

	while (<IN>) {
		if ($_ =~ m/^@/) {
			#################
			#################
			print OUT_SH $_;
			#################
			##################
			next;
		}
		my @line = split(/\s+/,$_);
		my $sam_flag = $line[1];
		if ($line[2] eq "*") {
			next;
		}
		if ( ($sam_flag == 99) || ($sam_flag == 355) || ($sam_flag == 147) || ($sam_flag == 403) || ($sam_flag == 83) || ($sam_flag == 339) || ($sam_flag == 163) || ($sam_flag == 419) ) {
			my ($read_id, $backsplice_candidate_name, $sam_start, $sam_cigar, $mate_sam_start) = ($line[0], $line[2], $line[3], $line[5], $line[7]);
			my ($pure_name, $e1, $e2);
			
			chomp($backsplice_candidate_name);

			if ($backsplice_candidate_name =~ m/^(.*)_E(\d+)B(\d+)$/) {
				$pure_name = $1;
				$e1 = $2;
				$e2 = $3;
			}
			my ($cigar_M) = $sam_cigar =~ m/(\d+)M/;
			my $transcript_name = $transcripts_seen{$pure_name};
			my ( $chr, $start, $E_lengths, $E_offs, $strand);
			if ($transcript_name =~ m/^(.*):(\d+)-\d+_.*_Exon_Lengths_(.*)_Offsets_(.*)_(\+|-)/) {
				$chr = $1;
				$start = $2;
				$E_lengths = $3;
				$E_offs = $4;
				$strand = $5;
			}
			
			#filter out reads mapped to the unexpected strand if strand specific mode turned on
			unless ($lib_type eq "fr-unstranded") {
				if ( ($sam_flag == 99) || ($sam_flag == 355) || ($sam_flag == 147) || ($sam_flag == 403) ) {
					unless ($strand eq $first_strand_expectation) {
						next;
					}
				} elsif ( ($sam_flag == 83) || ($sam_flag == 339) || ($sam_flag == 163) || ($sam_flag == 419) ) {
					unless ($strand eq $second_strand_expectation) {
						next;
					}
				}
			}
			
			my @E_lengths_ar = split(",",$E_lengths);

			#determine breakpoint position
			my $breakpoint_bp;
			if ($E_lengths_ar[$e1] < $READ_LENGTHS) {
				$breakpoint_bp = $E_lengths_ar[$e1];
			} else {
				$breakpoint_bp = $READ_LENGTHS-1;
			}
			
			#determine linear plus minus region length
			my $plus_minus_region_bp;
			if ($linear_ratio_info) {
				$plus_minus_region_bp = int ( ($breakpoint_bp / 2) + 0.5);
			}
			
			#test whether read crossed breakpoint
			if ( (($breakpoint_bp - $sam_start) >= ($MBSL - 1)) && (($sam_start + $cigar_M - $breakpoint_bp) >= ($MBSL - 1)) ) {
				my @E_offs_ar = split(",",$E_offs);
				my $purer_name;
				if ($pure_name =~ m/\d+-\d+_/) {
					#($purer_name) = $pure_name =~ m/chr.*:\d+-\d+_(.*;.*)_Exon_Lengths.*/;
					($purer_name) = $pure_name =~ m/^.*:\d+-\d+_(.*)_Exon_Lengths.*/;
				} else {
					$purer_name = $pure_name;
				}

				my $junction;
				if ($strand eq "+") {
					$junction = $purer_name . "_E" . ($e1+1) . "B" . ($e2+1);
				} elsif ($strand eq "-") {
					my $tot_exon_num = scalar(@E_lengths_ar);
					$junction = $purer_name . "_E" . ($tot_exon_num-$e2) . "B" . ($tot_exon_num-$e1);
				}
				
				my $genomic_start = $start + $E_offs_ar[$e2];
				my $genomic_end = $start + $E_offs_ar[$e1] + $E_lengths_ar[$e1];
								
				my $loci = "$chr\t$genomic_start\t$genomic_end";
				
				$name_reads{$backsplice_candidate_name}{$read_id}++;
				if (!exists $raw_to_output{$backsplice_candidate_name}) {
					$raw_to_output{$backsplice_candidate_name} = $junction;
				}
				if (!exists $raw_to_loci{$backsplice_candidate_name}) {
					$raw_to_loci{$backsplice_candidate_name} = $loci;
				}

				$loci_reads{$loci}{$read_id}{junc}++;
				$loci_tx_exons_names{$loci}{$junction}++;
				
				##extract and hash the beg/end segments
				my ($read_seq, $read_qual) = ($line[9], $line[10]);
				my $end_seg_length = $breakpoint_bp - $sam_start + 1;
				my $beg_seg_length = length($read_seq) - $end_seg_length;
				
				my ($end_seg_seq, $end_seg_qual, $beg_seg_seq, $beg_seg_qual);
				
				#if direction of mapped read is "left" (rev compl)
				if (($sam_flag == 83) || ($sam_flag == 147) || ($sam_flag == 339) || ($sam_flag == 403) ) {
					$read_seq = rev_compl($read_seq);
					$read_qual = rev_compl($read_qual);
				
					$beg_seg_seq = substr($read_seq, 0, $beg_seg_length);
					$beg_seg_qual = substr($read_qual, 0, $beg_seg_length);
					$end_seg_seq = substr($read_seq, $beg_seg_length, $end_seg_length);
					$end_seg_qual = substr($read_qual, $beg_seg_length, $end_seg_length);
				} else {
					$end_seg_seq = substr($read_seq, 0, $end_seg_length);
					$end_seg_qual = substr($read_qual, 0, $end_seg_length);
					$beg_seg_seq = substr($read_seq, $end_seg_length, $beg_seg_length);
					$beg_seg_qual = substr($read_qual, $end_seg_length, $beg_seg_length);
				}
				my $beg_line_2_4 = "$beg_seg_seq\n+\n$beg_seg_qual";
				my $end_line_2_4 = "$end_seg_seq\n+\n$end_seg_qual";
				
				my $in_pair_read_id = $in_pair{$sam_flag} . "&" . $read_id;
				my $backsplice_candidate_name_beg = $backsplice_candidate_name . "&beg";
				my $backsplice_candidate_name_end = $backsplice_candidate_name . "&end";
				
				$segment_fastq_seq{$beg_line_2_4}{$in_pair_read_id}{$backsplice_candidate_name_beg}++;
				$segment_fastq_seq{$end_line_2_4}{$in_pair_read_id}{$backsplice_candidate_name_end}++;

				$backsplice_candidate_name_breakpoint_bp{$backsplice_candidate_name} = $breakpoint_bp;
				
			#else if olf
			} elsif ( (scalar(@E_lengths_ar) == 1) || (($e1 == (scalar(@E_lengths_ar)-1)) && ($e2 == 0)) ){
				my @E_offs_ar = split(",",$E_offs);
				#if the read_id for this loci is already determined as junction, skip it
				my $genomic_start = $start + $E_offs_ar[$e2];
				my $genomic_end = $start + $E_offs_ar[$e1] + $E_lengths_ar[$e1];				
				my $loci = "$chr\t$genomic_start\t$genomic_end";
				if (exists $loci_reads{$loci}{$read_id}) {
					next;
				}
				
				my $direction;
				if (($sam_flag == 83) || ($sam_flag == 147) || ($sam_flag == 339) || ($sam_flag == 403)) {
					$direction = "left";
				} else {
					$direction = "right";
				}
				
				if ( ($direction eq "left") || (($linear_ratio_info) && ($direction eq "right")) ) {
					#check if read spans the first exon (crossed by BOUNDARY_THRESH) and does not exceed the junction breakpoint
					my $first_exon_end = $breakpoint_bp + $E_lengths_ar[0];
					if ( (($first_exon_end - $sam_start) >= ($MBSL - 1)) && (($sam_start - $breakpoint_bp) >= 1) ) {
						my $incr_string = $in_pair{$sam_flag} . "left";
						$potential_olf_reads{$backsplice_candidate_name}{$read_id}{$incr_string}++ unless ($direction eq "right");
						
						if ($linear_ratio_info) {
							$linear_loci_reads{$loci}{$read_id}{linF}++;
							$linear_name_reads{$backsplice_candidate_name}{$read_id}++;
						}

					}
				} 
				if ( ($direction eq "right") || (($linear_ratio_info) && ($direction eq "left")) ) {
					
					#dir is right; check if read spans the last exon (crosed by BOUNDARY_THRESH) OR mapped to annealed exon (cannot cross the breapoint)
					if (scalar(@E_lengths_ar) == 1) {
						if ($sam_start < ($breakpoint_bp - $MBSL)) {
							my $incr_string = $in_pair{$sam_flag} . "right";
							$potential_olf_reads{$backsplice_candidate_name}{$read_id}{$incr_string}++ unless ($direction eq "left");

							if ($linear_ratio_info) {
								$linear_loci_reads{$loci}{$read_id}{linF}++;
								$linear_name_reads{$backsplice_candidate_name}{$read_id}++;
							}
							
						}
						elsif ( ($sam_start > $breakpoint_bp) && (($sam_start - $mate_sam_start) >= $cigar_M) ) {
							my $incr_string = $in_pair{$sam_flag} . "right";
							$potential_olf_reads{$backsplice_candidate_name}{$read_id}{$incr_string}++ unless ($direction eq "left");
							
							if ($linear_ratio_info) {
								$linear_loci_reads{$loci}{$read_id}{linF}++;
								$linear_name_reads{$backsplice_candidate_name}{$read_id}++;
							}							
						}
					} elsif (scalar(@E_lengths_ar) > 1) {
						my $second_to_last_exon_end = $breakpoint_bp;
						for(my $i=0; $i < scalar(@E_lengths_ar) - 1; $i++) {
							$second_to_last_exon_end += $E_lengths_ar[$i];
						}
						my $sam_end = $sam_start + $cigar_M - 1;
						if ($sam_start < ($breakpoint_bp - $MBSL)) {
							my $incr_string = $in_pair{$sam_flag} . "right";
							$potential_olf_reads{$backsplice_candidate_name}{$read_id}{$incr_string}++ unless ($direction eq "left");
							
							if ($linear_ratio_info) {
								$linear_loci_reads{$loci}{$read_id}{linF}++;
								$linear_name_reads{$backsplice_candidate_name}{$read_id}++;
							}								
						}
						elsif ( (($sam_end - $second_to_last_exon_end) >= $MBSL) && ($sam_start > $breakpoint_bp) && (($sam_start - $mate_sam_start) >= $cigar_M) ) {
							my $incr_string = $in_pair{$sam_flag} . "right";
							$potential_olf_reads{$backsplice_candidate_name}{$read_id}{$incr_string}++ unless ($direction eq "left");
							
							if ($linear_ratio_info) {
								$linear_loci_reads{$loci}{$read_id}{linF}++;
								$linear_name_reads{$backsplice_candidate_name}{$read_id}++;
							}					
						}
					}
				}
			} elsif ($linear_ratio_info) {
				
				#test for reads falling in plus_minus_region
				if ( ( (($sam_start + $cigar_M) >= $plus_minus_region_bp) && (($sam_start + $cigar_M) <= $breakpoint_bp) ) || 
				     ( ($sam_start >= $breakpoint_bp) && ($sam_start <= ($breakpoint_bp + $plus_minus_region_bp)) )
					) {
				
					my @E_offs_ar = split(",",$E_offs);
					#if the read_id for this loci is already determined as junction, skip it
					my $genomic_start = $start + $E_offs_ar[$e2];
					my $genomic_end = $start + $E_offs_ar[$e1] + $E_lengths_ar[$e1];				
					my $loci = "$chr\t$genomic_start\t$genomic_end";
					
					$linear_name_reads{$backsplice_candidate_name}{$read_id}++;
					$linear_loci_reads{$loci}{$read_id}{linJ}++;
				
				}
			}
			print OUT_P $_;
		} else {
		}
	}
	close IN;

	close OUT_P;
	
	######
	close OUT_SH;
	#######
	
	my ($output_exon_ends_p1, $output_exon_ends_p2, $output_exon_ends_p3) = populate_and_extract_final_step_exon_ends_fasta(\%backsplice_candidate_name_breakpoint_bp, "bsj");
	my ($segments_fastq_file_p1, $segments_fastq_file_p2, $segments_fastq_file_p3) = print_segment_fastq_seq(\%segment_fastq_seq, "bsj");

	my %checked_junctions;
	print "Final round confirmation of backsplice segments...\n";
	map_and_examine_final_segments($output_exon_ends_p1, $segments_fastq_file_p1, $p_1_kmer, \%checked_junctions, "bsj");
	map_and_examine_final_segments($output_exon_ends_p2, $segments_fastq_file_p2, $p_2_kmer, \%checked_junctions, "bsj");	
	map_and_examine_final_segments($output_exon_ends_p3, $segments_fastq_file_p3, $p_3_kmer, \%checked_junctions, "bsj");	
	
	for my $backsplice_candidate_name (keys %potential_olf_reads) {
		for my $read_id (keys %{$potential_olf_reads{$backsplice_candidate_name}}) {
			my @read_info = (keys %{$potential_olf_reads{$backsplice_candidate_name}{$read_id}});
			next if (scalar(@read_info) != 2);
			@read_info = sort { substr($a,0,1) <=> substr($b,0,1) } @read_info;
			my ($key1, $key2) = ($read_info[0], $read_info[1]);
				
			my $comb1_right = "1right";
			my $comb1_left = "2left";
				
			my $comb2_right = "2right";
			my $comb2_left = "1left";
				
			my $found_pair = 0;
			if ( ($key1 eq $comb1_right) && ($key2 eq $comb1_left) ) {
				$found_pair = 1;
			} elsif ( ($key1 eq $comb2_left) && ($key2 eq $comb2_right) ) {
				$found_pair = 1;
			}
				
			if ($found_pair == 1) {
			
				my ($pure_name, $e1, $e2);

				if ($backsplice_candidate_name =~ m/^(.*)_E(\d+)B(\d+)$/) {
					$pure_name = $1;
					$e1 = $2;
					$e2 = $3;
				}

				my ($chr, $start, $E_lengths, $E_offs, $strand);
				my $transcript_name = $transcripts_seen{$pure_name};
				
				if ($transcript_name =~ m/^(.*):(\d+)-\d+_.*_Exon_Lengths_(.*)_Offsets_(.*)_(\+|-)/) {
					$chr = $1;
					$start = $2;
					$E_lengths = $3;
					$E_offs = $4;
					$strand = $5;
				}
			
				my @E_lengths_ar = split(",",$E_lengths);
				my @E_offs_ar = split(",",$E_offs);
				
				my $purer_name;
				if ($pure_name =~ m/\d+-\d+_/) {
					($purer_name) = $pure_name =~ m/^.*:\d+-\d+_(.*)_Exon_Lengths.*/;
				} else {
					$purer_name = $pure_name;
				}
				
				my $junction;
				if ($strand eq "+") {
					$junction = $purer_name . "_E" . ($e1+1) . "B" . ($e2+1);
				} elsif ($strand eq "-") {
					my $tot_exon_num = scalar(@E_lengths_ar);
					$junction = $purer_name . "_E" . ($tot_exon_num-$e2) . "B" . ($tot_exon_num-$e1);
				}
				
				my $genomic_start = $start + $E_offs_ar[$e2];
				my $genomic_end = $start + $E_offs_ar[$e1] + $E_lengths_ar[$e1];
								
				my $loci = "$chr\t$genomic_start\t$genomic_end";

				$name_reads{$backsplice_candidate_name}{$read_id}++;
				if (!exists $raw_to_output{$backsplice_candidate_name}) {
					$raw_to_output{$backsplice_candidate_name} = $junction;
				}
				if (!exists $raw_to_loci{$backsplice_candidate_name}) {
					$raw_to_loci{$backsplice_candidate_name} = $loci;
				}
				
				$loci_reads{$loci}{$read_id}{olf}++;
				$loci_tx_exons_names{$loci}{$junction}++;
			}
		}
	}

	my $filtered_proper_pair_out = $output_dir . "/candidate_circ_supporting_reads.sam";
	open OUT_PF, ">$filtered_proper_pair_out" or die "can't open $filtered_proper_pair_out\n";
	
	my $linear_out_temp;
	if ($linear_ratio_info) {
		$linear_out_temp = $output_dir . "/linear_supporting_reads_around_candidate_circ.temp";
		open OUT_LO, ">$linear_out_temp" or die "can't open $linear_out_temp\n";
	}
	
	my %filtered_name_reads;
	open IN, "<$proper_pair_out" or die "can't open $proper_pair_out\n";
	while (<IN>) {
		my @line = split("\t", $_);
		my ($read_id, $backsplice_candidate_name) = ($line[0], $line[2]);
		if (exists $name_reads{$backsplice_candidate_name}{$read_id}) {
			print OUT_PF $_;
			$filtered_name_reads{$backsplice_candidate_name}{$read_id}++;
		}
		
		if ($linear_ratio_info) {
			if (exists $linear_name_reads{$backsplice_candidate_name}{$read_id}) {
				if (exists $raw_to_loci{$backsplice_candidate_name}) {
					my $loci = $raw_to_loci{$backsplice_candidate_name};
					if (exists $loci_reads{$loci}{$read_id}{olf}) {
						delete $linear_loci_reads{$loci}{$read_id}{linF} if (exists $linear_loci_reads{$loci}{$read_id}{linF});
					} else {
						$line[2] = $raw_to_output{$backsplice_candidate_name};
						print OUT_LO join("\t", @line);
					}					
				}
			}
		}
	}
	close IN;
	close OUT_PF;
	if ($linear_ratio_info) {
		close OUT_LO;
	}
	
	#used by the -u/--unique-by option
	my %reads_loci;
	
	#Remove the unsatisfactory reads from filtered proper pair
	my $temp_filtered_proper_pair_out = $filtered_proper_pair_out . ".temp";
	system("mv $filtered_proper_pair_out $temp_filtered_proper_pair_out");
	open OUT_PF, ">$filtered_proper_pair_out" or die "can't open $filtered_proper_pair_out\n";
	open IN, "<$temp_filtered_proper_pair_out" or die "can't open $temp_filtered_proper_pair_out\n";
	while (<IN>) {
		my @line = split("\t", $_);
		my ($read_id, $backsplice_candidate_name) = ($line[0], $line[2]);
		my $loci = $raw_to_loci{$backsplice_candidate_name};
		
		unless ( ($filtered_name_reads{$backsplice_candidate_name}{$read_id} == $EXPECT_READIDS_PER_BSJ) && 
				((exists $loci_reads{$loci}{$read_id}{olf}) || 
				((exists $checked_junctions{$read_id}{$backsplice_candidate_name}{"1beg"}) && (exists $checked_junctions{$read_id}{$backsplice_candidate_name}{"1end"})) || 
				((exists $checked_junctions{$read_id}{$backsplice_candidate_name}{"2beg"}) && (exists $checked_junctions{$read_id}{$backsplice_candidate_name}{"2end"})) ) ) {
			
			delete $loci_reads{$loci}{$read_id} if (exists $loci_reads{$loci}{$read_id});
			next;
		}
		
		$reads_loci{$read_id}{$loci}++;
		
		$outputted_raw_names{$backsplice_candidate_name}++;

		print OUT_PF $_;
	}
	close IN;
	close OUT_PF;
	system("rm $temp_filtered_proper_pair_out");
	
	
	system("mv $filtered_proper_pair_out $temp_filtered_proper_pair_out");
	open OUT_PF, ">$filtered_proper_pair_out" or die "can't open $filtered_proper_pair_out\n";
	open IN, "<$temp_filtered_proper_pair_out" or die "can't open $temp_filtered_proper_pair_out\n";
	while (<IN>) {
		my @line = split("\t", $_);
		my ($read_id, $backsplice_candidate_name) = ($line[0], $line[2]);
		my $loci = $raw_to_loci{$backsplice_candidate_name};
	
		if ((scalar(keys %{$reads_loci{$read_id}})) > $unique_by) {
			unless (exists $loci_reads{$loci}{$read_id}{olf}) {
				delete $loci_reads{$loci}{$read_id} if (exists $loci_reads{$loci}{$read_id});
				delete $outputted_raw_names{$backsplice_candidate_name} if (exists $outputted_raw_names{$backsplice_candidate_name});
				next;
			}
		}
		
		$line[2] = $raw_to_output{$backsplice_candidate_name};
		print OUT_PF join("\t", @line);
	}	
	
	close OUT_PF;
	close IN;
	system("rm $temp_filtered_proper_pair_out");
	
	

	my $filtered_sam_headers = $final_remapping . "/filtered_sam_headers";
	open OUT_FSH, ">$filtered_sam_headers" or die "can't open $filtered_sam_headers\n";
	open IN_SH, "<$sam_headers" or die "can't open $sam_headers\n";
	while (<IN_SH>) {
		my @line = split("\t", $_);
		my $target;
		if ($line[1] =~ /SN:(.*)/) {
			$target = $1;
			if (exists $outputted_raw_names{$target}) {
				$line[1] = "SN:" . $raw_to_output{$target};
				print OUT_FSH join("\t", @line);
			}
		}
	}
	close IN_SH;
	close OUT_FSH;
	
	system("mv $filtered_proper_pair_out $temp_filtered_proper_pair_out");
	system("cat $filtered_sam_headers $temp_filtered_proper_pair_out > $filtered_proper_pair_out");
	system("rm $temp_filtered_proper_pair_out");
	
	
	my $output_bsj_fa = $output_dir . "/candidate_circ_junctions.fa";
	open OUT_BFA, ">$output_bsj_fa" or die "can't open $output_bsj_fa\n";
	open IN_TFA, "<$target_potential_backsplice_fasta" or die "can't open $target_potential_backsplice_fasta\n";
	while (<IN_TFA>) {
		chomp;
		if ($_ =~ /^>(.*)/) {
			my $target = $1;
			if (exists $outputted_raw_names{$target}) {
				my $seq = <IN_TFA>;
				print OUT_BFA ">$raw_to_output{$target}\n$seq";
			}
		}
	}
	close OUT_BFA;
	close IN_TFA;
	
	my $linear_out;
	if ($linear_ratio_info) {
		$linear_out = $output_dir . "/linear_supporting_reads_around_candidate_circ.sam";
		system("mv $linear_out_temp $linear_out");
		system("rm $linear_out_temp");
	}	

	my %unique_reads;
	my %highest_readcount;
	
	my %loci_linear;
	
	my $output_candidate_circ_genomic_junctions = "$output_dir/candidate_circ_junctions.bed";
	my $num_circRNAs=0;
	for my $loci (keys %loci_reads) {
		if (%{$loci_reads{$loci}}) {
			my @reads = keys (%{$loci_reads{$loci}});
			
			my @tx_exon_names = sort keys (%{$loci_tx_exons_names{$loci}});		
			my $tx_exon_names_var = join("=", @tx_exon_names);
			
			#find how many junc reads (if any) and olf reads (if any) there are
			my ($junc_reads,$olf_reads) = (0,0);
			for my $read_id (@reads) {
				if (exists $loci_reads{$loci}{$read_id}{junc}) {
					$junc_reads++;
					$unique_reads{$read_id}{$tx_exon_names_var}++;
				} elsif (exists $loci_reads{$loci}{$read_id}{olf}) {
					$olf_reads++;
					$unique_reads{$read_id}{$tx_exon_names_var}++;
				}
			}
			my ($linear_J_reads, $linear_F_reads) = (0,0);
			if ($linear_ratio_info) {
				if (exists $linear_loci_reads{$loci}) {
					for my $lin_read_id (keys %{$linear_loci_reads{$loci}}) {
						if (exists $linear_loci_reads{$loci}{$lin_read_id}{linJ}) {
							$linear_J_reads++;
						} elsif (exists $linear_loci_reads{$loci}{$lin_read_id}{linF}) {
							$linear_F_reads++;
						}
					}
				}
			}			
						
			my $tot_reads = 0;
			$tot_reads = $junc_reads + $olf_reads;
			
			my $score;
			if ($olf_reads == 0) {
				$score = 1;
			} elsif ($junc_reads == 0) {
				$score = 0;
			} else {
				$score = $junc_reads / $tot_reads;
			}

			my $linear_reads;
			if ($score > 0.5) {
				$linear_reads = $linear_J_reads;
			} else {
				$linear_reads = $linear_F_reads;
			}
			
			my $final_name;
			if ($linear_ratio_info) {
				$final_name = $tot_reads . ":" . $tx_exon_names_var . ":" . $linear_reads;
			} else {
				$final_name =  $tot_reads . ":" . $tx_exon_names_var;
			}
			if ( ($tot_reads != 0) || ($linear_ratio_info) ) {
				$loci_linear{$loci}{circ} = $tot_reads;
				$loci_linear{$loci}{lin} = $linear_reads;
				$num_circRNAs++;
				$highest_readcount{"$loci\t$final_name\t$score"} = $tot_reads;
			}
			
		}
	}
		
	my $output_and_raw_junction_names;
	if ($linear_ratio_info) {
		$output_and_raw_junction_names = $output_dir . "/linear_info.tsv";
	} else {
		$output_and_raw_junction_names = $output_dir . "/junction_names_output_to_raw.tsv";
	}
	open OUT_ORN, ">$output_and_raw_junction_names" or die "can't open $output_and_raw_junction_names\n";
	
	if ($linear_ratio_info) {
		print OUT_ORN "output_name\t" . "raw_name\t" . "genomic_loci\t" . "linear_reads\t" . "approx_linear_ratio";
		if ($equiv_fsj) {
			print OUT_ORN "\t" . "nearest_fsj_equiv_to_bsj\n";
		} else {
			print OUT_ORN "\n";
		}
	} else {
		print OUT_ORN "output_name\t" . "raw_name\t" . "genomic_loci";
		if ($equiv_fsj) {
			print OUT_ORN "\t" . "nearest_fsj_equiv_to_bsj\n";
		} else {
			print OUT_ORN "\n";
		}
	}

	for my $backsplice_candidate_name (keys %raw_to_output) {
		if (!exists $raw_to_loci{$backsplice_candidate_name}) {
			print "$backsplice_candidate_name does not exist in raw_to_loci\n";
			next;
		}
		my $loci = $raw_to_loci{$backsplice_candidate_name};
		my @loci_ar = split("\t", $loci);
		my $loci_final = $loci_ar[0] . ":" . $loci_ar[1] . "-" . $loci_ar[2];
		if ($linear_ratio_info) {
			if (exists $loci_linear{$loci}) {
				my $circ_linear = $loci_linear{$loci}{circ} + $loci_linear{$loci}{lin};
				my $approx_linear_ratio;
				if ($circ_linear == 0) {
					$approx_linear_ratio = 0;
				} else {
					$approx_linear_ratio = $loci_linear{$loci}{circ} / $circ_linear;
				}
				print OUT_ORN $raw_to_output{$backsplice_candidate_name} . "\t" . $backsplice_candidate_name . "\t" .
							$loci_final . "\t" . $loci_linear{$loci}{lin} . "\t" . $approx_linear_ratio;
				if ($equiv_fsj) {
					my $nearest_fsj = ".";
				
					my ($pure_name, $bsj);
					if ($backsplice_candidate_name =~ m/^(.*)_(E\d+B\d+)$/) {
						$pure_name = $1;
						$bsj = $2;
					}
					if (exists $bsjs_with_equiv_fsj{$pure_name}) {
						if (exists $bsjs_with_equiv_fsj{$pure_name}{$bsj}) {
							$nearest_fsj = $bsjs_with_equiv_fsj{$pure_name}{$bsj};
						}
					}
					print OUT_ORN "\t" . $nearest_fsj . "\n";
				} else {
					print OUT_ORN "\n";
				}
			}
		} else {
			print OUT_ORN $raw_to_output{$backsplice_candidate_name} . "\t" . $backsplice_candidate_name . "\t" .
							$loci_final;
							
			if ($equiv_fsj) {
				my $nearest_fsj = ".";
			
				my ($pure_name, $bsj);
				if ($backsplice_candidate_name =~ m/^(.*)_(E\d+B\d+)$/) {
					$pure_name = $1;
					$bsj = $2;
				}
				if (exists $bsjs_with_equiv_fsj{$pure_name}) {
					if (exists $bsjs_with_equiv_fsj{$pure_name}{$bsj}) {
						$nearest_fsj = $bsjs_with_equiv_fsj{$pure_name}{$bsj};
					}
				}
				print OUT_ORN "\t" . $nearest_fsj . "\n";
			} else {
				print OUT_ORN "\n";
			}
		}
	}
	close OUT_ORN;

	open OUT_H, ">$output_candidate_circ_genomic_junctions" or die "can't open $output_candidate_circ_genomic_junctions\n";
	for my $line (sort { $highest_readcount{$b} <=> $highest_readcount{$a} } keys %highest_readcount) {
		print OUT_H $line . "\n";
	}
	close OUT_H;
	
	print "Number of backspliced unique reads: " . scalar(keys %unique_reads) . "\n";
	
	#my ($NBPF_related, $not_NBPF_related) = (0,0);
	#for my $read (keys %unique_reads) {
	#	my $related_flag=0;
	#	for my $final_name (keys %{$unique_reads{$read}}) {
	#		if ($final_name =~ m/NBPF/) {
	#			$NBPF_related++;
	#			$related_flag=1;
	#			last;
	#		}
	#	}
	#	if ($related_flag==0) {
	#		$not_NBPF_related++;
	#	}
	#}
	#print "Name contains the word NBPF: $NBPF_related\n";
	#print "does not contain the word NBPF: $not_NBPF_related\n";
	print "\nNumber of candidate circular RNA junctions: $num_circRNAs\n";
	
	
	if ($linear_ratio_info) {
		system("rm $proper_pair_out");
	}
	
	my $outputted_raw_names_list = $mapping_dir . "/outputted_raw_names.lst";
	open OUT, ">$outputted_raw_names_list" or die "can't open $outputted_raw_names_list\n";

	for my $backsplice_candidate_name (keys %outputted_raw_names) {
		print OUT $backsplice_candidate_name . "\t" . $outputted_raw_names{$backsplice_candidate_name} . "\n";
	}
	close OUT;
	
	if (!$full_length) {
		return;
	}

	for my $backsplice_candidate_name (keys %outputted_raw_names) {
		
		my ($pure_name, $bsj);
		if ($backsplice_candidate_name =~ m/^(.*)_(E\d+B\d+)$/) {
			$pure_name = $1;
			$bsj = $2;
		}
		my ($e, $b);
		if ($bsj =~ m/E(\d+)B(\d+)/) {
			$e = $1;
			$b = $2;
		}
		
		my $transcript_name = $transcripts_seen{$pure_name};
		my $strand;
		if ($transcript_name =~ m/^.*:\d+-\d+_.*_Exon_Lengths_.*_Offsets_.*_(\+|-)/) {
			$strand = $1;
		}
		
		$raw_name_strand{$backsplice_candidate_name} = $strand;

		if ( ($e - $b) == 0 ) {
			next;
		}

		$transcripts_bsjs{$transcript_name}{$bsj}++;	
	}
}


sub final_mapping_and_output_fsjs {
	my $full_length_isoforms_fa_index = $target_potential_fli_fasta . ".index";
	system("$kallisto index -k $final_kallisto_k -i $full_length_isoforms_fa_index $target_potential_fli_fasta");
	
	my $final_remapping = $full_length_dir;
	
	my $proper_pair_out = $final_remapping . "/fsj_reads_mapped_proper_pair.sam";

	my $sam_headers = $final_remapping . "/sam_headers";
	open OUT_SH, ">$sam_headers" or die "can't open $sam_headers\n";
	
	open OUT_P, ">$proper_pair_out" or die "can't open $proper_pair_out\n";
	
	if ( ($low_num) && ($low_num_sam == 1)) {
		open IN, "<$low_num_temp_sam" or die "can't open $low_num_temp_sam\n";
	} elsif ($low_num) {	
		open IN, "$kallisto pseudo -i $full_length_isoforms_fa_index -o $final_remapping -t 1 --pseudobam $potential_fsj_reads_1 $potential_fsj_reads_2 |"  or die "can't run kallisto pseudo\n";
	} else {
		open IN, "$kallisto pseudo -i $full_length_isoforms_fa_index -o $final_remapping -t $threads_kallisto --pseudobam $potential_fsj_reads_1 $potential_fsj_reads_2 |"  or die "can't run kallisto pseudo\n";
	}

	
	my %segment_fastq_seq;
	
	#used for extracting final step exon_ends_fasta
	my %full_length_isoform_name_breakpoint_bp;
	
	my ($first_strand_expectation, $second_strand_expectation);
	if ($lib_type eq "fr-firststrand") {
		$first_strand_expectation = "-";
		$second_strand_expectation = "+";
	} elsif ($lib_type eq "fr-secondstrand") {
		$first_strand_expectation = "+";
		$second_strand_expectation = "-";
	}

	my %checked_junctions;
	
	while (<IN>) {
		if ($_ =~ m/^@/) {
			print OUT_SH $_;
			next;
		}
		my @line = split(/\s+/,$_);
		next if ($line[2] eq "*");
		my $sam_flag = $line[1];
		if ( ($sam_flag == 99) || ($sam_flag == 355) || ($sam_flag == 147) || ($sam_flag == 403) || ($sam_flag == 83) || ($sam_flag == 339) || ($sam_flag == 163) || ($sam_flag == 419) ) {
		my ($read_id, $full_length_isoform_name, $sam_start, $sam_cigar, $mate_sam_start) = ($line[0], $line[2], $line[3], $line[5], $line[7]);
			
			my $direction;
			if (($sam_flag == 83) || ($sam_flag == 147) || ($sam_flag == 339) || ($sam_flag == 403)) {
				$direction = "left";
			} else {
				$direction = "right";
			}
			
			#must be inward facing
			if ($direction eq "right") {
				next unless ($mate_sam_start >= $sam_start);
			} else {
				next unless ($sam_start >= $mate_sam_start);
			}
			
			print OUT_P $_;

			my ($cigar_M) = $sam_cigar =~ m/(\d+)M/;
			
			my $strand = $fli_strands{$full_length_isoform_name};
			
			#filter out reads mapped to the unexpected strand if strand specific mode turned on
			unless ($lib_type eq "fr-unstranded") {
				if ( ($sam_flag == 99) || ($sam_flag == 355) || ($sam_flag == 147) || ($sam_flag == 403) ) {
					unless ($strand eq $first_strand_expectation) {
						next;
					}
				} elsif ( ($sam_flag == 83) || ($sam_flag == 339) || ($sam_flag == 163) || ($sam_flag == 419) ) {
					unless ($strand eq $second_strand_expectation) {
						next;
					}
				}
			}

			for my $fsj (keys %{$fli_breakpoints{$full_length_isoform_name}}) {
				my $breakpoint_bp = $fli_breakpoints{$full_length_isoform_name}{$fsj};
				#test whether read crossed breakpoint
				if ( (($breakpoint_bp - $sam_start) >= ($MBSL - 1)) && (($sam_start + $cigar_M - $breakpoint_bp) >= ($MBSL - 1)) ) {
				
					my ($E, $B) = $fsj =~ m/E(\d+)B(\d+)/;
					my $isoform_fsj = $full_length_isoform_name . "_" . $fsj;
					
					
						##extract and hash the beg/end segments
						my ($read_seq, $read_qual) = ($line[9], $line[10]);
						my $end_seg_length = $breakpoint_bp - $sam_start + 1;
						my $beg_seg_length = length($read_seq) - $end_seg_length;
				
						my ($end_seg_seq, $end_seg_qual, $beg_seg_seq, $beg_seg_qual);

						if ($direction eq "left") {

							$read_seq = rev_compl($read_seq);
							$read_qual = rev_compl($read_qual);
				
							$beg_seg_seq = substr($read_seq, 0, $beg_seg_length);
							$beg_seg_qual = substr($read_qual, 0, $beg_seg_length);
							$end_seg_seq = substr($read_seq, $beg_seg_length, $end_seg_length);
							$end_seg_qual = substr($read_qual, $beg_seg_length, $end_seg_length);
						} else {
							$end_seg_seq = substr($read_seq, 0, $end_seg_length);
							$end_seg_qual = substr($read_qual, 0, $end_seg_length);
							$beg_seg_seq = substr($read_seq, $end_seg_length, $beg_seg_length);
							$beg_seg_qual = substr($read_qual, $end_seg_length, $beg_seg_length);
						}
						my $beg_line_2_4 = "$beg_seg_seq\n+\n$beg_seg_qual";
						my $end_line_2_4 = "$end_seg_seq\n+\n$end_seg_qual";
				
						my $in_pair_read_id = $in_pair{$sam_flag} . "&" . $read_id;
						
						my $full_length_isoform_name_beg = $isoform_fsj . "&beg";
						my $full_length_isoform_name_end = $isoform_fsj . "&end";
				
						$segment_fastq_seq{$beg_line_2_4}{$in_pair_read_id}{$full_length_isoform_name_beg}++;
						$segment_fastq_seq{$end_line_2_4}{$in_pair_read_id}{$full_length_isoform_name_end}++;

						$full_length_isoform_name_breakpoint_bp{$full_length_isoform_name}{$fsj} = $breakpoint_bp;
				}
			} 
		}
	}
	close IN;

	close OUT_P;

	my ($output_exon_ends_p1, $output_exon_ends_p2, $output_exon_ends_p3) = populate_and_extract_final_step_exon_ends_fasta(\%full_length_isoform_name_breakpoint_bp, "fsj");
	my ($segments_fastq_file_p1, $segments_fastq_file_p2, $segments_fastq_file_p3) = print_segment_fastq_seq(\%segment_fastq_seq, "fsj");

	print "Final round confirmation of fsj segments...\n";
	map_and_examine_final_segments($output_exon_ends_p1, $segments_fastq_file_p1, $p_1_kmer, \%checked_junctions, "fsj");
	map_and_examine_final_segments($output_exon_ends_p2, $segments_fastq_file_p2, $p_2_kmer, \%checked_junctions, "fsj");	
	map_and_examine_final_segments($output_exon_ends_p3, $segments_fastq_file_p3, $p_3_kmer, \%checked_junctions, "fsj");	
	
	open INV, ">>$invalidated_flis" or die "can't open $invalidated_flis\n";
	
	my %confirmed_alt_fsj_reads;
	my %confirmed_alt_fsjs;
	for my $full_length_isoform_name (keys %fli_breakpoints) {
		my @empty_fsjs;
		for my $fsj (keys %{$fli_breakpoints{$full_length_isoform_name}}) {
			my $isoform_fsj = $full_length_isoform_name . "_" . $fsj;
			if (exists $checked_junctions{$isoform_fsj}){
					for my $read_id (keys %{$checked_junctions{$isoform_fsj}}) {
						if ( ((exists $checked_junctions{$isoform_fsj}{$read_id}{"1beg"}) && (exists $checked_junctions{$isoform_fsj}{$read_id}{"1end"})) || 
						((exists $checked_junctions{$isoform_fsj}{$read_id}{"2beg"}) && (exists $checked_junctions{$isoform_fsj}{$read_id}{"2end"})) ) {
							$confirmed_alt_fsj_reads{$read_id}++;
						} else {
							delete $checked_junctions{$isoform_fsj}{$read_id};
						}
					}
				if (scalar(keys %{$checked_junctions{$isoform_fsj}}) == 0) {
					push @empty_fsjs, $fsj;
				}
			} else {
				push @empty_fsjs, $fsj;
			}
		}
		if (scalar(@empty_fsjs) > 0) {
			print INV "$full_length_isoform_name invalidated due to no paired reads mapped at " . join(", ", @empty_fsjs) . "\n";
		} else {
			$confirmed_alt_fsjs{$full_length_isoform_name}++;
		}
	}
	
	close INV;
	
	#get raw_name to loci and raw_name to output sequence
	my %raw_name_out;
	my %og_to_raw;
	my $output_and_raw_junction_names = $output_dir . "/junction_names_output_to_raw.tsv";
	open IN, "<$output_and_raw_junction_names" or die "can't open $output_and_raw_junction_names\n";
	my $header = <IN>;
	while (<IN>) {
		chomp;
		my @line = split("\t", $_);
		my ($output_name,$raw_name,$genomic_loci) = ($line[0],$line[1],$line[2]);
		$raw_name_out{$raw_name}{o} = $output_name;
		$raw_name_out{$raw_name}{g} = $genomic_loci;
		$og_to_raw{"$output_name\t$genomic_loci"} = $raw_name;
	}
	close IN;
	
	my %raw_to_output_name;
	#get the junction to outward ratios of raw names
	my $output_candidate_circ_genomic_junctions = "$output_dir/candidate_circ_junctions.bed";
	open IN, "<$output_candidate_circ_genomic_junctions" or die "can't open $output_candidate_circ_genomic_junctions\n";
	while (<IN>) {
		chomp;
		my @line = split("\t", $_);
		my ($chr,$start,$end,$info,$ratio) = @line;
		my $loci = "$chr:$start-$end";
		my ($output_name) = $info =~ m/^\d+:(.*)$/;
		my @output_name_ar = split("=", $output_name);
		for my $indiv_output_name (@output_name_ar) {
			my $raw_name = $og_to_raw{"$indiv_output_name\t$loci"};
			if (!$raw_name) {
				print "Cannot locate $indiv_output_name in tsv\n";
			}
			$raw_to_output_name{$raw_name} = $indiv_output_name;
			$raw_name_out{$raw_name}{r} = $ratio;
		}
	}
	close IN;	
	
	#output full length isoform fasta
	my $fli_output_fasta = $full_length_dir . "/full_length_isoforms_unmerged.fa";
	open OUT, ">$fli_output_fasta" or die "can't open $fli_output_fasta\n";
	#open IN, "<$target_potential_backsplice_fasta" or die "can't open $target_potential_backsplice_fasta\n";
	open IN, "<$target_potential_backsplice_fasta_no_copied_seq" or die "can't open $target_potential_backsplice_fasta_no_copied_seq\n";
	while (<IN>) {
		chomp;
		if ($_ =~ m/^>(.*)/) {
			next unless (exists $outputted_raw_names{$1});
			next if (! exists $raw_name_out{$1}{r});
			my $raw_name = $1;
			my ($E,$B);
			if ($raw_name =~ m/.*_E(\d+)B(\d+)/) {
				$E = $1;
				$B = $2;
			}
			
			
			my $seq = <IN>;
			chomp($seq);
			
			#my $assumed_breakpoint = $READ_LENGTHS - 1;
			#if (length($seq) <= $assumed_breakpoint) {
			#	my $corrected_seq = substr($seq, (length($seq)/2));
			#	print OUT ">" . $raw_to_output_name{$raw_name} . "\tC\n$corrected_seq\n";			
			#	next;
			#}
			#my $minus_annealed = substr($seq, $assumed_breakpoint);
			#my $length_minus_annealed = length($minus_annealed);
			#if ($length_minus_annealed >= $assumed_breakpoint) {
			#	print OUT ">" . $raw_to_output_name{$raw_name} ."\tC\n$minus_annealed\n";
			#} else {
			#	my $corrected_seq = substr($seq, (length($seq)/2));
			#	print OUT ">" . $raw_to_output_name{$raw_name} . "\tC\n$corrected_seq\n";
			#}
			
			print OUT ">" . $raw_to_output_name{$raw_name} . "\tC\n$seq\n";
		}
	}
	close IN;
	close OUT;
	
	my %final_fli_output;
	my %total_fli_loci;
	my %total_fli;
	
	
	##append originnal annotated linear transcript isoforms
	open IN, "<$fa" or die "can't open $fa\n";
	open OUT, ">>$fli_output_fasta" or die "can't open $fli_output_fasta\n";
	while (<IN>) {
		chomp;
		if ($_ =~ m/^>(.*:\d+-\d+)_(.*)_Exon_Lengths_.*_Offsets_.*_(\+|-)/) {
			my $seq = <IN>;
			#chomp($seq);
			
			my $transcript_genomic_loci = $1;
			my $pure_name = $2;
			my $strand = $3;

			my ($transcript_name) = $_ =~ m/^>(.*)/;
			
			my $raw_name;
			if (exists $transcripts_seen{$transcript_name}) {
				$raw_name = $transcript_name;
			} else {
				$raw_name = $pure_name;
			}
			
			print OUT ">$raw_name\n$seq";
			
			
			
			my $final_output_string = $raw_name . "\t" . $raw_name . "\t" . $transcript_genomic_loci 
				. "\t" . "." . "\t" . "." . "\t" . $strand
				. "\t.\t.\n";

			$final_fli_output{$final_output_string} = 0;

			$total_fli_loci{lin}{$transcript_genomic_loci}++;
			$total_fli{lin}++;
			
		}
	}
	close IN;
	close OUT;
	
	open INV, ">>$invalidated_flis" or die "can't open $invalidated_flis\n";
	
	for my $raw_name (keys %outputted_raw_names) {
		my ($E,$B);
		if ($raw_name =~ m/.*_E(\d+)B(\d+)/) {
			$E = $1;
			$B = $2;
		}
		
		my $bsj_read_count = $outputted_raw_names{$raw_name}/2;
		
		if (! exists $raw_name_out{$raw_name}{r}) {
			print INV $raw_name_out{$raw_name}{o} . "\t" . $raw_name . "\t" . $raw_name_out{$raw_name}{g} . "\t" 
					. "have inconsistent output loci and output backsplice transcript\n";
			next;
		}
		my $final_output_string = $raw_name_out{$raw_name}{o} . "\t" . $raw_name . "\t" . $raw_name_out{$raw_name}{g} 
			. "\t" . $bsj_read_count . "\t" . $raw_name_out{$raw_name}{r} . "\t" . $raw_name_strand{$raw_name} 
			. "\t.\t.\n";
		$final_fli_output{$final_output_string} = $bsj_read_count;
		
		$total_fli_loci{circ}{$raw_name_out{$raw_name}{g}}++;
		$total_fli{circ}++;
	}

	my %confirmed_alt_flis;
	
	my %alt_raw_to_output_name;
	for my $full_length_isoform_name (keys %fli_breakpoints) {

		next unless (exists $confirmed_alt_fsjs{$full_length_isoform_name});
		
		my ($raw_name, $pure_name);
		if ($full_length_isoform_name =~ m/(.*)(_E\d+B\d+)_\d+A\d+.*$/) { #modified from .*_\d+A\d$ to _\d+A\d+.*$
			$pure_name = $1;
			$raw_name = $1 . $2;
		} elsif ($full_length_isoform_name =~ m/(.*?)_\d+A\d+.*$/) { #change from .*\d+A\d+$ to \d+A\d+.*$
			$pure_name = $1;
			$raw_name = $1;
		}

		my $transcript_name = $transcripts_seen{$pure_name};
		my ($transcript_genomic_loci, $E_lengths) = $transcript_name =~ m/^(.*:\d+-\d+)_.*_Exon_Lengths_(.*)_Offsets_.*_/;
		my @E_lengths_ar = split(",", $E_lengths);
		my $num_exons = scalar(@E_lengths_ar);		
		
		my $strand = $fli_strands{$full_length_isoform_name};
		#my %final_output_adj_fsjs;
		my %final_output_alt_fsjs;
		for my $fsj (keys %{$fli_breakpoints{$full_length_isoform_name}}) {
			my $isoform_fsj = $full_length_isoform_name . "_" . $fsj;
			my $fsj_read_count = scalar(keys %{$checked_junctions{$isoform_fsj}});
			
			my ($fsj_E, $fsj_B) = $fsj =~ m/E(\d+)B(\d+)/;
			my ($final_fsj_E, $final_fsj_B);
			if ($strand eq "+") {
				$final_fsj_E = $fsj_E+1;
				$final_fsj_B = $fsj_B+1;
			} else {
				$final_fsj_E = $num_exons-$fsj_B;
				$final_fsj_B = $num_exons-$fsj_E;
			}
			
			my $final_adj_or_alt_fsj = $final_fsj_E . "_" . $final_fsj_B;
			$final_output_alt_fsjs{$final_adj_or_alt_fsj} = $fsj_read_count;
		}

		
		my $alt_info;
		if ($full_length_isoform_name =~ m/.*_E\d+B\d+_(\d+A\d+.*)$/) { #move .* to the end
			$alt_info = $1;
		} elsif ($full_length_isoform_name =~ m/.*?_(\d+A\d+.*)$/) { #move .* to the end 
			$alt_info = $1;
		}
		
		
		my @alt_info_ar;
		my @final_alt_info_ar;
		my $alt_out_name;
		if ($alt_info) {
			@alt_info_ar = split("_", $alt_info);

			for my $alt (@alt_info_ar) {
				my ($alt_E, $alt_B);
				if ($alt =~ m/(\d+)A(\d+)/) {
					$alt_E = $1;
					$alt_B = $2;
				}
				my ($final_alt_E, $final_alt_B);
				if ($strand eq "+") {
					$final_alt_E = $alt_E+1;
					$final_alt_B = $alt_B+1;
				} else {
					$final_alt_E = $num_exons-$alt_B;
					$final_alt_B = $num_exons-$alt_E;
				}
				my $final_alt = $final_alt_E . "A" . $final_alt_B;
				push @final_alt_info_ar, $final_alt;
			}
			@final_alt_info_ar = sort { ($a =~ /(\d+)/)[0] <=> ($b =~ /(\d+)/)[0] } @final_alt_info_ar;
			
			my $out_raw_name;
			if ($raw_name =~ m/_E\d+B\d+$/) {
				$out_raw_name = $raw_name_out{$raw_name}{o};
			} else {
				$out_raw_name = $raw_name;
			}
			
			$alt_out_name = $out_raw_name . "_" . join("_",@final_alt_info_ar);
			
			$alt_raw_to_output_name{$full_length_isoform_name} = $alt_out_name;
		}
		
		my $bsj_read_count;
		if ($raw_name =~ m/_E\d+B\d+$/) {
			$bsj_read_count = $outputted_raw_names{$raw_name}/2;
		} else {
			$bsj_read_count = ".";
		}
		
		my $final_output_string;
		if ($alt_info) {
			$final_output_string .= $alt_out_name;
		} else {
			$final_output_string .= $raw_name_out{$raw_name}{o};
		}
		if ($raw_name =~ m/_E\d+B\d+$/) {
			if (! exists $raw_name_out{$raw_name}{r}) {
				print INV $raw_name_out{$raw_name}{o} . "\t" . $raw_name . "\t" . $raw_name_out{$raw_name}{g} . "\t" 
						. "have inconsistent output loci and output backsplice transcript\n";
						next;
			}
		}
		
		my ($output_genomic_loci, $output_ratio);
		if ($raw_name =~ m/_E\d+B\d+$/) {
			$output_genomic_loci = $raw_name_out{$raw_name}{g};
			$output_ratio = $raw_name_out{$raw_name}{r};
			
			$total_fli_loci{circ}{$output_genomic_loci}++;
			$total_fli{circAlt}++;
		} else {
			$output_genomic_loci = $transcript_genomic_loci;
			$output_ratio = ".";
			
			$total_fli_loci{lin}{$output_genomic_loci}++;
			$total_fli{linAlt}++;
		}
		
		$final_output_string .= "\t" . $full_length_isoform_name . "\t" . $output_genomic_loci
				. "\t" . $bsj_read_count . "\t" . $output_ratio . "\t" . $strand 
				. "\t";
		
		
				$confirmed_alt_flis{$full_length_isoform_name}++;

		my @final_alt_fsjs;
		my @final_alt_fsjs_rcs;

		if ($alt_info) {
			for my $final_fsj (sort { ($a =~ /(\d+)/)[0] <=> ($b =~ /(\d+)/)[0] } keys %final_output_alt_fsjs) {
				push @final_alt_fsjs, $final_fsj;
				push @final_alt_fsjs_rcs, $final_output_alt_fsjs{$final_fsj};
			}
			$final_output_string .= join("|",@final_alt_fsjs) . "\t" . join("|",@final_alt_fsjs_rcs) . "\n";


		} else {
			$final_output_string .= ".\t.\n";
		
		}
		
		my $sort_val;
		if ($bsj_read_count eq ".") {
			$sort_val = 0;
		} else {
			$sort_val = $bsj_read_count;
		}
		$final_fli_output{$final_output_string} = $sort_val;

	}
	
	open IN, "<$proper_pair_out" or die "can't open $proper_pair_out\n";
	my $alt_fsj_supporting_reads_temp = $full_length_dir . "/temp_alt_fsj_supporting_reads.sam";
	my $alt_fsj_supporting_reads = $full_length_dir . "/full_length_isoforms_alt_fsj_supporting_reads_unmerged.sam";
	open OUT, ">$alt_fsj_supporting_reads_temp" or die "can't open $alt_fsj_supporting_reads_temp\n";
	while (<IN>) {
		my @line = split("\t", $_);
		if (exists $confirmed_alt_fsj_reads{$line[0]}) {
			my $full_length_isoform_name = $line[2];
			#next if (exists $invalidated_flis{$full_length_isoform_name});
			next unless (exists $confirmed_alt_flis{$full_length_isoform_name});
			$line[2] = $alt_raw_to_output_name{$full_length_isoform_name};
			print OUT join("\t", @line);
		}
	}
	close IN;
	close OUT;
	
	system("rm $proper_pair_out");

	my $filtered_sam_headers = $final_remapping . "/filtered_sam_headers";
	open OUT_FSH, ">$filtered_sam_headers" or die "can't open $filtered_sam_headers\n";
	open IN_SH, "<$sam_headers" or die "can't open $sam_headers\n";
	while (<IN_SH>) {
		my @line = split("\t", $_);
		chomp(@line);
		my ($target, $length);
		if ($line[1] =~ /SN:(.*)$/) {
			$target = $1;
			if ($line[2] =~ /LN:(\d+)$/) {
				$length = $1;
			}
			if (exists $confirmed_alt_flis{$target}) {
				my $out_name = $alt_raw_to_output_name{$target};
				
				print OUT_FSH "\@SQ\t" . "SN:$out_name\t" . "LN:$length\n";
			}
		}
	}
	close IN_SH;
	close OUT_FSH;
	
	system("cat $filtered_sam_headers $alt_fsj_supporting_reads_temp > $alt_fsj_supporting_reads");
	system("rm $alt_fsj_supporting_reads_temp");
	
	
	open OUT, ">>$fli_output_fasta" or die "can't open $fli_output_fasta\n";
	open IN, "<$target_potential_fli_fasta" or die "can't open $target_potential_fli_fasta\n";
	while (<IN>) {
		chomp;
		if ($_ =~ m/^>(.*)/) {
			my $full_length_isoform_name = $1;
			next unless (exists $confirmed_alt_flis{$full_length_isoform_name});
			
			if ($_ =~ m/_E\d+B\d+/) {
				print OUT ">" . $alt_raw_to_output_name{$full_length_isoform_name} . "\tC\n";
			} else {
				print OUT ">" . $alt_raw_to_output_name{$full_length_isoform_name} . "\n";
			}
			
			my $seq = <IN>;
			print OUT $seq;
		}
	}
	close OUT;
	close IN;
	
	
	my $fli_output = $full_length_dir . "/full_length_isoforms_unmerged.tsv";
	open OUT, ">$fli_output" or die "can't open $fli_output\n";
	print OUT "output_name\t" . "raw_name\t" . "genomic_loci\t" . "back_splice_junction_read_count\t" . "cross_bsj_to_olf_reads_ratio\t" 
				. "strand\t" . "alternative_forward_splice_junctions\t" . "alternative_forward_splice_junctions_read_count\n";
	for my $line (sort { $final_fli_output{$b} <=> $final_fli_output{$a} } keys %final_fli_output) {
		print OUT $line;
	}
	close OUT;
	close INV;
	
	my $total_genomic_loci=0;
	for my $type (keys %total_fli_loci) {
		$total_genomic_loci += scalar(keys %{$total_fli_loci{$type}});
	}
	my $total_fli=0;
	for my $type (keys %total_fli) {
		$total_fli += $total_fli{$type};
	}
	
	print "Total number of genomic loci: $total_genomic_loci\n";
	for my $type (qw(lin circ)) {
		my $name;
		if ($type eq "lin") {
			$name = "linear (input annotation): ";
		} elsif ($type eq "circ") {
			$name = "circular (back-splicing junction): ";
		}
	
		print "\t" . $name . scalar(keys %{$total_fli_loci{$type}}) . "\n";
	}
	
	
	print "\n" . "Number of full-length isoforms: $total_fli\n";
	for my $type (qw(lin circ linAlt circAlt)) {
		my $name;
		if ($type eq "lin") {
			$name = "linear (input annotation): ";
		} elsif ($type eq "circ") {
			$name = "circular (back-splicing junction): ";
		} elsif ($type eq "linAlt") {
			$name = "linear alternative (exon skipping): ";
		} elsif ($type eq "circAlt") {
			$name = "circular alternative (exon skipping): ";
		}
		my $val;
		if (!exists $total_fli{$type}) {
			$val=0;
		} else {
			$val = $total_fli{$type};
		}
		print "\t" . $name . $val . "\n";
	}
	
	merge_full_length_output($fli_output_fasta, $alt_fsj_supporting_reads, $fli_output);
}


sub merge_full_length_output {
	my ($fli_output_fasta_unmerged, $alt_fsj_supporting_reads_unmerged, $fli_output_unmerged) = @_;
	
	my %isoform_to_strand;
	open my $IN, "<$fli_output_unmerged" or die "can't open $fli_output_unmerged\n";
	my $header = <$IN>;
	while (<$IN>) {
		my @line = split("\t", $_);
		my ($isoform, $strand) = ($line[0], $line[5]);
		$isoform_to_strand{$isoform} = $strand;
	}
	close $IN;
	
	my %unique_fli_seqs;
	open $IN, "<$fli_output_fasta_unmerged" or die "can't open $fli_output_fasta_unmerged\n";
	while (<$IN>) {
		chomp;
		if ($_ =~ /^>(\S+)/) {
			my $fli_name = $1;
			my $type;
			if ($_ =~ /\s+C$/) {
				$type = "C";
			} else {
				$type = "L";
			}
			my $seq = <$IN>;
			chomp($seq);
			$unique_fli_seqs{$type}{$isoform_to_strand{$fli_name}}{$seq}{$fli_name}++;
		}
	}
	close $IN;
	
	my %unmerged_to_merged_name;
	
	for my $type (keys %unique_fli_seqs) {
		for my $strand (keys %{$unique_fli_seqs{$type}}) {
			for my $seq (keys %{$unique_fli_seqs{$type}{$strand}}) {
				my @merged_fli_names = sort keys %{$unique_fli_seqs{$type}{$strand}{$seq}};
				my $merged_fli_name = join(",", @merged_fli_names);
				for my $fli_name (@merged_fli_names) {
					$unmerged_to_merged_name{$fli_name} = $merged_fli_name;
				}
			}
		}
	}
	
	output_merged_fli_fasta($fli_output_fasta_unmerged, \%unmerged_to_merged_name);
	output_merged_fli_alt_fsj_sam($alt_fsj_supporting_reads_unmerged, \%unmerged_to_merged_name);
	output_merged_fli_output($fli_output_unmerged, \%unmerged_to_merged_name);
}

sub output_merged_fli_output {
	my ($fli_output_unmerged, $unmerged_to_merged_name) = @_;
	
	my $fli_output = $output_dir . "/full_length_isoforms.tsv";
	open my $OUT, ">$fli_output" or die "can't open $fli_output\n";
	open my $IN, "<$fli_output_unmerged" or die "can't open $fli_output_unmerged\n";
	
	my %already_outputted;
	
	my $header = <$IN>;
	my @header_ar = split("\t", $header);
	splice @header_ar, 1, 1;
	print $OUT join("\t", @header_ar);
	
	my %merged_fli_strands;
	
	while (<$IN>) {
		my @line = split("\t", $_);
		
		my $strand = $line[5];
		
		splice @line, 1, 1;
		
		my $fli_name = $line[0];
		my $merged_name = $unmerged_to_merged_name->{$fli_name};
		
		if (exists $already_outputted{$merged_name}) {
			next;
		}
		
		$line[0] = $merged_name;
		print $OUT join("\t", @line);
		
		$merged_fli_strands{$merged_name} = $strand;
		
		$already_outputted{$merged_name}++;
	}
	
	close $OUT;
	close $IN;
	
	my $fli_output_fasta_ref_plus_strand = $full_length_dir . "/full_length_isoforms.RefPlusStrand.fa";
	open $IN, "<$fli_output_fasta_ref_plus_strand" or die "can't open $fli_output_fasta_ref_plus_strand\n";
	my $fli_output_fasta = $output_dir . "/full_length_isoforms.fa";
	open $OUT, ">$fli_output_fasta" or die "can't open $fli_output_fasta\n";
	while (<$IN>) {
		if ($_ =~ /^>(\S+)/) {
			my $merged_name = $1;
			my $header = $_;
			my $strand = $merged_fli_strands{$merged_name};
			
			my $seq = <$IN>;
			chomp($seq);
			if ($strand eq "-") {
				$seq = rev_compl($seq);
			}
			
			print $OUT $header . "$seq\n";		
		}
	}
	close $OUT;
	close $IN;
}

sub output_merged_fli_alt_fsj_sam {
	my ($alt_fsj_supporting_reads_unmerged, $unmerged_to_merged_name) = @_;
	
	my $alt_fsj_supporting_reads = $output_dir . "/full_length_isoforms_alt_fsj_supporting_reads.sam";
	open my $OUT, ">$alt_fsj_supporting_reads" or die "can't open $alt_fsj_supporting_reads\n";
	open my $IN, "<$alt_fsj_supporting_reads_unmerged" or die "can't open $alt_fsj_supporting_reads_unmerged\n";
	my %already_outputted_header;
	my %chosen_unmerged;
	while (<$IN>) {
		chomp;
		if ($_ =~ /^@/) {
			if ($_ =~ /SN:(\S+).*LN:(\d+)/) {
				my $fli_name = $1;
				my $length = $2;
				
				my $merged_name = $unmerged_to_merged_name->{$fli_name};
				if (!$merged_name) {
					print "No merged name for $fli_name\n";
					next;
				}
				if (exists $already_outputted_header{$merged_name}) {
					next;
				}
				
				print $OUT "\@SQ\t" . "SN:$merged_name\tLN:$length\n";
				
				$already_outputted_header{$merged_name}++;
				$chosen_unmerged{$fli_name}++;
				
				next;
			}
		}
		my @line = split("\t", $_);
		
		if (exists $chosen_unmerged{$line[2]}) {
			my $chosen_fli_name = $line[2];
			$line[2] = $unmerged_to_merged_name->{$chosen_fli_name};
			print $OUT join("\t", @line) . "\n";
		} else {
			next;
		}
	}
	close $IN;
	close $OUT;
	
}

sub output_merged_fli_fasta {
	my ($fli_output_fasta_unmerged, $unmerged_to_merged_name) = @_;
	

	my $fli_output_fasta_ref_plus_strand = $full_length_dir . "/full_length_isoforms.RefPlusStrand.fa";
	open my $OUT_R, ">$fli_output_fasta_ref_plus_strand" or die "can't open $fli_output_fasta_ref_plus_strand\n";

	#my $fli_output_fasta = $output_dir . "/full_length_isoforms.fa";
	#open my $OUT, ">$fli_output_fasta" or die "can't open $fli_output_fasta\n";
	open my $IN, "<$fli_output_fasta_unmerged" or die "can't open $fli_output_fasta_unmerged\n";
	
	my %already_outputted;
	while (<$IN>) {
		if ($_ =~ /^>(\S+)/) {
			my $fli_name = $1;
			
			my $strand = $fli_strands{$fli_name};
			
			my $merged_name = $unmerged_to_merged_name->{$fli_name};
			
			if (exists $already_outputted{$merged_name}) {
				next;
			}
			
			if ($_ =~ /\s+C$/) {
				print $OUT_R ">" . $merged_name . "\t" . "C\n";

			} else {
				print $OUT_R ">" . $merged_name . "\n";

			}
			my $seq = <$IN>;
			print $OUT_R $seq;
			
			$already_outputted{$merged_name}++;
			
		}
	}
	close $IN;
	close $OUT_R;
}


sub map_reads_to_targets {

	#get the hash for existing transcripts first
	open IN, "<$fa" or die "can't open $fa\n";
	while (<IN>) {
		chomp;
		if ($_ =~ m/^.*:\d+-\d+_(.*)_Exon_Lengths_.*_Offsets_.*_/) {
			my $pure_name = $1;
			my ($transcript_name) = $_ =~ m/^>(.*)/;

			if ( (exists $transcripts_seen{$pure_name}) && ($transcripts_seen{$pure_name} ne $transcript_name) ) {
				$transcripts_seen{$transcript_name} = $transcript_name;
				$pure_name = $transcript_name;
			}
			else {
				$transcripts_seen{$pure_name} = $transcript_name;
			}
		}
		
	}
	close IN;

	open IN, "$kallisto pseudo -i $all_possible_bsj_targets_fa_index -o $mapping_dir --single -t $threads_kallisto -l $kallisto_l -s 20 --pseudobam $fq1 $fq2 |" or die "can't run kallisto pseudo\n";

	my %transcripts;
	my %read_ids;
	my %read_lengths;
	
	print "Mapping reads to all_possible bsj targets... ";
	system("date");
	
	while (<IN>) {
		next if (/^@/);
		my @line = split("\t", $_);
		next if ($line[1] == 4);
		
		$read_ids{$line[0]}++;

		if ($line[2] =~ m/(.*)_(E\d+B\d+)$/) {
			$transcripts{$1}{$2}++;
		}
		$read_lengths{length($line[10])}++;

	}
	close IN;

		
	print "Finished mapping. ";
	system("date");

	#Invalidate transcripts with > max_mappable_bsjs
	open OUT, ">$mapping_dir/invalidated_too_many_mappable_bsjs" or die "can't open $mapping_dir/invalidated_too_many_mappable_bsjs\n";
	print OUT "Alter the --max-mappable-bsjs option to re-include the invalidated transcripts with too many mappable bsjs below:\n";
	for my $raw_name (keys %transcripts) {
		my $info = $transcripts_seen{$raw_name};
		my $num_mappable_bsjs = scalar(keys %{$transcripts{$raw_name}});
		if ($num_mappable_bsjs >= $max_mappable_bsjs) {
			delete $transcripts{$raw_name};
			print OUT "$raw_name invalidated because it has $num_mappable_bsjs (more than $max_mappable_bsjs" . ")" . " mappable bsjs. ";
			print OUT "Entry: $info\n\n";
		}
	}
	close OUT;
	
	my $most_common_matched_length;
	for my $read_length (sort { $read_lengths{$b} <=> $read_lengths{$a} } keys %read_lengths) {
		$most_common_matched_length = $read_length;
		last;
	}
	if (!$most_common_matched_length) {
		$most_common_matched_length = 101;
	}
	
	
	my $junc_reads = scalar(keys %read_ids);
	my $low_junc_reads;
	if ($junc_reads <= $low_num_reads) {
		$low_junc_reads++;
	}
	
	open IN, "$kallisto pseudo -i $FirstandLastExons_entities_fa_index -o $mapping_dir -t $threads_kallisto --pseudobam $fq1 $fq2 | " or die "can't open kallisto pseudo\n";
	
	print "Mapping reads to FirstandLastExon entities... ";
	system("date");


	
	while (<IN>) {
		next if (/^@/);
		my @line = split("\t",$_);
		if ( ($line[1] == 99) || ($line[1] == 355) || ($line[1] == 163) || ($line[1] == 419) ) {			
			

		    if ( ( ($line[3] - $line[7]) >= $most_common_matched_length) ) {
		    	$read_ids{$line[0]}++;

				if ($line[2] =~ m/(.*)_(E\d+B\d+)$/) {
					$transcripts{$1}{$2}++;
				}
				
				$read_lengths{length($line[10])}++;

		    } elsif ($low_junc_reads) {
				$low_num_read_ids{$line[0]}++;
			}
		} elsif ( ($line[1] == 83) || ($line[1] == 339) || ($line[1] == 147) || ($line[1] == 403) ) {
			
			if ( ( ($line[7] - $line[3]) >= $most_common_matched_length) ) {
		    	$read_ids{$line[0]}++;

				if ($line[2] =~ m/(.*)_(E\d+B\d+)$/) {
					$transcripts{$1}{$2}++;
				}
				$read_lengths{length($line[10])}++;
			} elsif ($low_junc_reads) {
				$low_num_read_ids{$line[0]}++;
			}
		}
		
	}
	close IN;
	#close OUT;
	$num_reads = scalar(keys %read_ids);
	#print "Num_reads: " . $num_reads . "\n";
	#browse(\%read_ids);
	$num_transcripts = scalar(keys %transcripts);
	#print "Num_tx: " . $num_transcripts . "\n";
	#browse(\%transcripts);
	if ( ($num_reads <= $low_num_reads) && ($num_transcripts <= $low_num_transcripts) ) {
		$low_num++;
		
		for my $read_id (keys %low_num_read_ids) {
			$read_ids{$read_id}++;
		}
	}
	
	print "Finished mapping. ";
	system("date");
	
	for my $read_length (sort { $read_lengths{$b} <=> $read_lengths{$a} } keys %read_lengths) {
		$most_common_matched_length = $read_length;
		last;
	}
	if (!$most_common_matched_length) {
		$most_common_matched_length = 101;
	}
	
	print "Extracting possible bsj reads... ";
	system("date");
	extract_reads($potential_backsplice_reads_1, $fq1, $potential_backsplice_reads_2, $fq2, \%read_ids);
	print "Finished extracting. ";
	system("date");

	print "Extracting potential backsplice transcripts... ";
	system("date");
	extract_backsplice_transcripts(\%transcripts, $most_common_matched_length);
	print "Finished extracting. ";
	system("date");
	
	if ($equiv_fsj) {
		print "Collecting equiv fsjs... ";
		system("date");
		collect_equiv_fsjs(\%transcripts);
		print "Finished collecting. ";
		system("date");
	}
	
	###
	if ($low_num) {
		my $backsplice_transcripts_fa_index = $target_potential_backsplice_fasta . ".index";
		system("$kallisto index -k $final_kallisto_k -i $backsplice_transcripts_fa_index $target_potential_backsplice_fasta");
		system("$kallisto pseudo -i $backsplice_transcripts_fa_index -o $mapping_dir -t $threads_kallisto --pseudobam $potential_backsplice_reads_1 $potential_backsplice_reads_2 > $low_num_temp_sam");
		check_temp_sam($low_num_temp_sam);
	}
	###
	
	return $most_common_matched_length;
}

sub check_temp_sam {
	my $sam_file = shift;
	open my $IN, "<$sam_file"or die "can't open $sam_file\n";
	while (<$IN>) {
		if ($_ =~ /^@/) {
			next;
		}
		my @line = split("\t", $_);
		if ($line[1]) {
			if ($line[1] =~ /^\d+$/) {
				$low_num_sam = 1;
				last;
			}
		}
	}
	close $IN;
}

sub map_reads_to_targets_fsj {
	my $fsj_targets_fa_index = $fsj_targets_fa . ".index";

	system("$kallisto index -i $fsj_targets_fa_index -k $kallisto_k $fsj_targets_fa");

	open IN, "$kallisto pseudo -i $fsj_targets_fa_index -o $full_length_dir --single -t $threads_kallisto -l $kallisto_l -s 20 --pseudobam $fq1 $fq2 |" or die "can't run kallisto pseudo\n";
	
	my %fsjs;
	my %read_ids;
	my %read_lengths;

	print "Mapping reads to fsj targets... ";
	system("date");
	
	while (<IN>) {
		next if (/^@/);
		my @line = split("\t", $_);
		next if ($line[1] == 4);
		
		$read_ids{$line[0]}++;

		if ($line[2] =~ m/(.*)_(E\d+B\d+)$/) {
			$fsjs{$1}{$2}++;
			
		}
		$read_lengths{length($line[10])}++;

	}
	close IN;

	$num_reads = scalar(keys %read_ids);
	$num_transcripts = scalar(keys %fsjs);
	if ( ($num_reads <= $low_num_reads) && ($num_transcripts <= $low_num_transcripts) ) {
		$low_num++;
		for my $read_id (keys %low_num_read_ids) {
			$read_ids{$read_id}++;
		}
	}
	
	print "Finished mapping. ";
	system("date");

	print "Extracting potential full-length circRNA isoforms... ";
	system("date");
	extract_full_length_isoforms(\%fsjs);
	print "Finished extracting. ";
	system("date");
	
	
	print "Extracting potential fsj reads... ";
	system("date");
	extract_reads($potential_fsj_reads_1, $fq1, $potential_fsj_reads_2, $fq2, \%read_ids);
	print "Finished extracting. ";
	system("date");
	
	if ($low_num) {
		my $full_length_isoforms_fa_index = $target_potential_fli_fasta . ".index";
		system("$kallisto index -k $final_kallisto_k -i $full_length_isoforms_fa_index $target_potential_fli_fasta");
		system("$kallisto pseudo -i $full_length_isoforms_fa_index -o $full_length_dir -t $threads_kallisto --pseudobam $potential_fsj_reads_1 $potential_fsj_reads_2 > $low_num_temp_sam");
		check_temp_sam($low_num_temp_sam);
	}

}

sub extract_full_length_isoforms {
	my $fsjs = shift;
	
	open INV, ">$invalidated_flis" or die "can't open $invalidated_flis\n";
	
	my %full_length_to_extract;
	
	#%$fsjs contains ALL mappable alt fsjs for each transcript, including those outside of bsj
	for my $raw_name (keys %$fsjs) {
		my $transcript_name = $transcripts_seen{$raw_name};
		my $E_lengths;
		my $strand;
		if ($transcript_name =~ m/^.*:\d+-\d+_.*_Exon_Lengths_(.*)_Offsets_.*_(\+|-)/) {
			$E_lengths = $1;
			$strand = $2;
		}
		my @E_lengths_ar = split(",",$E_lengths);
		my $l_B = 0;
		my $l_E = scalar(@E_lengths_ar)-1;

		#define E and B exons, and also if "b" (back-splice) or "l" (linear original annotated)
		gather(\$raw_name, $fsjs, $l_E, $l_B, "l", \$transcript_name, $strand, \%full_length_to_extract);
		for my $bsj (keys %{$transcripts_bsjs{$transcript_name}}) {

			
			my ($E, $B);
			if ($bsj =~ m/E(\d+)B(\d+)/) {
				$E = $1;
				$B = $2;
			}
			gather(\$raw_name, $fsjs, $E, $B, "b", \$transcript_name, $strand, \%full_length_to_extract);
		}
	}

	extract(\%full_length_to_extract);
	
	close INV;

}

sub gather {
	my ($raw_name, $fsjs, $E, $B, $l_or_b, $transcript_name, $strand, $full_length_to_extract) = @_;

	my $raw_name_bsj;
	if ($l_or_b eq "l") {
		$raw_name_bsj = $$raw_name;
	} elsif ($l_or_b eq "b") {
		$raw_name_bsj = $$raw_name . "_E$E" . "B$B";
	}
	
	$full_length_to_extract->{$$transcript_name}->{$raw_name_bsj}++;	
	$fli_strands{$raw_name_bsj} = $strand;
	
	my %alt_fsjs;
	my $num_alt_fsjs = 0;
	for my $E_and_B (keys %{$fsjs->{$$raw_name}}) {
		my ($fsj_E, $fsj_B) = $E_and_B =~ m/E(\d+)B(\d+)/;

		if ( (($fsj_B >= $B) && ($fsj_B <= $E)) && (($fsj_E >= $B) && ($fsj_E <= $E)) ) {
			$alt_fsjs{$fsj_E}{$fsj_B}++;
			$num_alt_fsjs++;
		}
	}
	
	return if ($num_alt_fsjs == 0);
	
	if ($num_alt_fsjs > $max_mappable_alt_fsjs) {
		print INV "$raw_name_bsj invalidated due to more than $max_mappable_alt_fsjs ($num_alt_fsjs) number of alt fsj\n";
		return;
	}
	
	my %set_alt_fsjs;
	my %considered;
	my $alt_fsj_in_question;
	while (scalar(keys %considered) != $num_alt_fsjs) {
		OUTER:for my $fsj_E (sort {$a <=> $b} keys %alt_fsjs) {
			for my $fsj_B (sort {$a <=> $b} keys %{$alt_fsjs{$fsj_E}}) {
				my $E_and_B = "E$fsj_E" . "B$fsj_B";
				if (!exists $considered{$E_and_B}) {
					$alt_fsj_in_question = $E_and_B;
					last OUTER;
				}
			}
		}
		for my $fsj_E (sort {$a <=> $b} keys %alt_fsjs) {
			for my $fsj_B (sort {$a <=> $b} keys %{$alt_fsjs{$fsj_E}}) {
				my $E_and_B = "E$fsj_E" . "B$fsj_B";
				if ($E_and_B eq $alt_fsj_in_question) {
					push @{$set_alt_fsjs{$alt_fsj_in_question}}, $E_and_B;
					$considered{$alt_fsj_in_question}++;
				} elsif ( alt_fsj_coexist(\$E_and_B, \$alt_fsj_in_question) ) {
					push @{$set_alt_fsjs{$alt_fsj_in_question}}, $E_and_B;
					$considered{$E_and_B}++;
				}
			}
		}
	}

	my %considered_chains;
	my %combo_tested;
	for my $set (keys %set_alt_fsjs) {
		my @set = @{$set_alt_fsjs{$set}};

		
		my %set_considered;
		while (scalar(keys %set_considered) != scalar(@set)) {
			for (my $i=0; $i < @set; $i++) {
				
				
				my @chain;
			
				for (my $j=0; $j < $i; $j++) {
					if (!@chain) {
						push @chain, $set[$j];
						$set_considered{$set[$j]}++;
					} elsif (alt_fsj_coexist(\$set[$j], \$chain[-1]) ) {
						push @chain, $set[$j];
						$set_considered{$set[$j]}++;
					}
				}
				
				
				push @chain, $set[$i];
				$set_considered{$set[$i]}++;
				
				for (my $k=$i+1; $k < @set; $k++) {
					if (alt_fsj_coexist(\$set[$k], \$chain[-1]) ) {
						push @chain, $set[$k];
						$set_considered{$set[$k]}++;
					}
				}
				
				my $chain_str = join(" ", @chain);
			
				
				if (exists $considered_chains{$chain_str}) {

					next;
				}

				my @all_combo = findcombo(\@chain);
				for my $combo (@all_combo) {
					next if (exists $combo_tested{$combo});
					
					my @combo_ar = split(' ', $combo);
					$combo_tested{$combo}++;
					
					if (scalar(@combo_ar) >= $max_alt_fsj_per_isoform) {
						print INV "$combo for $raw_name_bsj invalidated due to having more than $max_alt_fsj_per_isoform connected alt fsjs\n";
						next;
					}
					
					my $inval;
					my @As;
					my $prev_B=0;
					for my $alt_E_and_B (@combo_ar) {
						my ($alt_E, $alt_B) = $alt_E_and_B =~ m/E(\d+)B(\d+)/;
						
						if ( ($alt_E - $prev_B) < 0) {
							$inval++;
						}
						
						my $A = $alt_E . "A" . $alt_B;
						push @As, $A;
						
						$prev_B = $alt_B;
					}
					
					if ($inval) {
						print INV "$combo for $raw_name_bsj invalidated due to alt fsjs that cannot coexist\n";
						next;
					}
					

					my $alt_name_bsj = $raw_name_bsj . "_" . join("_",@As);
					$full_length_to_extract->{$$transcript_name}->{$alt_name_bsj}++;
					$fli_strands{$alt_name_bsj} = $strand;
				}
			
				$considered_chains{$chain_str}++;
			}
		}
	}
}

sub extract {
	my $full_length_to_extract = shift;

	open IN, "<$fa" or die "can't open $fa\n";
	my %alt_fli_seq;
	while (<IN>) {
		chomp;
		if ($_ =~ /^>(.*)/) {
			my $transcript_name = $1;
			
			if (!exists $full_length_to_extract->{$transcript_name}) {
				next;
			}
			
			my $E_lengths;
			if ($transcript_name =~ m/^.*:\d+-\d+_.*_Exon_Lengths_(.*)_Offsets_.*_/) {
					$E_lengths = $1;
			}
			my @E_lengths_ar = split(",",$E_lengths);

			my $seq = <IN>;
			chomp($seq);
			
			for my $fli_name (keys %{$full_length_to_extract->{$transcript_name}}) {
				#my ($raw_name, $bsj) = $fli_name =~ m/(.*)_(E\d+B\d+)/;		
				
				my $alt_info;
				
				my ($E, $B);
				if ($fli_name =~ m/.*_E(\d+)B(\d+)/) { #Warning what if there exists gene of name E*B*
					$E = $1;
					$B = $2;
					
					if ($fli_name =~ m/.*_E\d+B\d+_(\d+A\d+.*)$/) { #move .* to the back 
						$alt_info = $1;
					}
					
				} else {
					$E = scalar(@E_lengths_ar)-1;
					$B = 0;
					
					if ($fli_name =~ m/.*?_(\d+A\d+(_\d+A\d+)*)$/) { #modified regex to match mutliple *A* 
						$alt_info = $1;
					}
				}
				
				my ($e1, $e2) = ($E, $B);
				
				my @alt_info_ar;
				if ($alt_info) {
					#@alt_info_ar = split(",", $alt_info);
					@alt_info_ar = split("_", $alt_info);
				}
				
				#if (!@alt_info_ar) {
					#print "No alt_info_ar for $fli_name\n";
					#browse(\@alt_info_ar);
					#}
				
				my $running_total = 0;
				my $isoform_sequence;
				if ($e2 != 0) {
					for (my $i=0; $i < $e2; $i++) {
						$running_total += $E_lengths_ar[$i];
					}
				}
				
				my @fli_exon_nums;
				my @fli_exon_lengths;
				
				for (my $i = $e2; $i <= $e1; $i++) {

					my $alt_exception=0;
					for my $alt_E_and_B (@alt_info_ar) {
						my ($alt_E, $alt_B) = $alt_E_and_B =~ m/(\d+)A(\d+)/;
						if ( ($i > $alt_E) && ($i < $alt_B) ) {
							$alt_exception=1;
							last;
						}					
					}
					if ($alt_exception==1) {
						$running_total += $E_lengths_ar[$i];
						next;
					}
					
					push @fli_exon_nums, $i;
					push @fli_exon_lengths, $E_lengths_ar[$i];
					
					$isoform_sequence .= substr($seq,$running_total,$E_lengths_ar[$i]);
					$running_total += $E_lengths_ar[$i];
					
				}
				
				if ($alt_info) {
					$alt_fli_seq{$fli_name} = $isoform_sequence;
				} else {
					#print OUT ">$fli_name\n$isoform_sequence\n";
				}
				
				my $tot_length;
				for (my $i=0; $i < scalar(@fli_exon_nums)-1; $i++) {
					$tot_length += $fli_exon_lengths[$i];
					
					#exception1
					if ( ($fli_exon_lengths[$i] < $MBSL) || ($fli_exon_lengths[$i+1] < $MBSL) ) {
						#print "Exception 1 $fli_name at $i or " . ($i+1) . "\n";
						next;
					}
				
					#exception2
					if ( ($fli_exon_lengths[$i] + $fli_exon_lengths[$i+1]) < ($kallisto_k + 1) ) {
						#print "Exception 2 $fli_name at $i or " . ($i+1) . "\n";
						next;
					}
					
					#non-alternative
					if ( ($fli_exon_nums[$i+1] - $fli_exon_nums[$i]) == 1) {
						next;
					}
					
					$fli_breakpoints{$fli_name}{"E$fli_exon_nums[$i]" . "B$fli_exon_nums[$i+1]"} = $tot_length;
				}
			}

							
			
		}
	}
	close IN;

	my %eliminated_alt_flis;
	
	check_identical_to_outputted_bst(\%alt_fli_seq, \%eliminated_alt_flis);
	check_identical_to_orig_transcript(\%alt_fli_seq, \%eliminated_alt_flis);
	
	
	open OUT, ">$target_potential_fli_fasta" or die "can't open $target_potential_fli_fasta\n";

	for my $fli_name (keys %alt_fli_seq) {
		if (exists $eliminated_alt_flis{$fli_name}) {
			delete $fli_breakpoints{$fli_name};
			print INV $fli_name . " invalidated due to having the exact same sequence as the following fli(s) with all adj fsjs:" . $eliminated_alt_flis{$fli_name} . "\n";
		} else {
			print OUT ">$fli_name\n$alt_fli_seq{$fli_name}\n";
		}
	}
	close OUT;
}

sub check_identical_to_orig_transcript {
	my ($alt_fli_seq, $eliminated_alt_flis) = @_;
	open IN, "<$fa" or die "can't open $fa\n";
	
	my %ot;
	while (<IN>) {
		chomp;
		if ($_ =~ m/^.*:\d+-\d+_(.*)_Exon_Lengths_.*_Offsets_.*_/) {
			my $seq = <IN>;
			
			my $pure_name = $1;
			
			my ($transcript_name) = $_ =~ m/^>(.*)/;
			
			my $raw_name;
			if (exists $transcripts_seen{$transcript_name}) {
				$raw_name = $transcript_name;
			} else {
				$raw_name = $pure_name;
			}

			push @{$ot{$seq}}, $raw_name;
			
		}
	}
	close IN;
	
	for my $fli_name (keys %{$alt_fli_seq}) {
		if (exists $ot{$alt_fli_seq->{$fli_name}}) {
			$eliminated_alt_flis->{$fli_name} .= " ";
			$eliminated_alt_flis->{$fli_name} .= join(" ", @{$ot{$alt_fli_seq->{$fli_name}}});
		}
	} 		
}

sub check_identical_to_outputted_bst {
	my ($alt_fli_seq, $eliminated_alt_flis) = @_;
	
	my %bst;
	
	#open IN, "<$target_potential_backsplice_fasta" or die "can't open $target_potential_backsplice_fasta\n";
	open IN, "<$target_potential_backsplice_fasta_no_copied_seq" or die "can't open $target_potential_backsplice_fasta_no_copied_seq\n";
	while (<IN>) {
		chomp;
		if ($_ =~ m/^>(.*)/) {
			next unless (exists $outputted_raw_names{$1});

			my $raw_name = $1;

			my $seq = <IN>;
			chomp($seq);
			
			#my $assumed_breakpoint = $READ_LENGTHS - 1;
			#if (length($seq) <= $assumed_breakpoint) {
			#	my $corrected_seq = substr($seq, (length($seq)/2));		
			#	push @{$bst{$corrected_seq}}, $raw_name;
			#	next;
			#}
			#my $minus_annealed = substr($seq, $assumed_breakpoint);
			#my $length_minus_annealed = length($minus_annealed);
			#if ($length_minus_annealed >= $assumed_breakpoint) {
			#	push @{$bst{$minus_annealed}}, $raw_name;
			#} else {
			#	my $corrected_seq = substr($seq, (length($seq)/2));
			#	push @{$bst{$corrected_seq}}, $raw_name;
			#}
			
			push @{$bst{$seq}}, $raw_name;
		
		}
	}
	close IN;
	
	for my $fli_name (keys %{$alt_fli_seq}) {
		if (exists $bst{$alt_fli_seq->{$fli_name}}) {
			$eliminated_alt_flis->{$fli_name} .= " "; 
			$eliminated_alt_flis->{$fli_name} .= join(" ", @{$bst{$alt_fli_seq->{$fli_name}}});
		}
	} 
}	
			
sub alt_fsj_coexist {
	my ($alt_fsj_1, $alt_fsj_2) = @_;
	
	my ($E_1, $B_1) = $$alt_fsj_1 =~ m/E(\d+)B(\d+)/;
	my ($E_2, $B_2) = $$alt_fsj_2 =~ m/E(\d+)B(\d+)/;
	
	my $coexist;
	if ( ($E_1 <= $E_2) && ($B_1 <= $E_2) ) {
		$coexist = 1;
	} elsif ( ($E_2 <= $E_1) && ($B_2 <= $E_1) ) {
		$coexist = 1;
	}
	if ($coexist) {
		return 1;
	} else {
		return 0;
	}
}


sub findcombo {

   my ($list) = @_;
   my (@combo, $str, $i, $j);

   my $size = @{$list};
   
   my @final_combo;

   for ($i = 0; $i < 2**$size; $i++) {
      $str = sprintf("%*.*b", $size, $size, $i);
      @combo = ();
      for ($j = 0; $j < $size; $j++) {
         if (substr($str, $j, 1)) { 
         	push (@combo, $list->[$j]);
         }
      }
      if (@combo) {
     	 push @final_combo, join(' ', @combo);
     }
   }
   return @final_combo;
}

sub collect_equiv_fsjs {
	my $transcripts = shift;
	my $each_side_length = $kallisto_k - $MBSL;
	open IN, "<$fa" or die "can't open $fa\n";
	while (<IN>) {
		chomp;
		if ($_ =~ /^>(.*)/) {

			my $transcript_name = $1;
			my ($E_lengths, $pure_name);
			if ($transcript_name =~ m/^.*:\d+-\d+_(.*)_Exon_Lengths_(.*)_Offsets_.*_/) {
					$pure_name = $1;
					$E_lengths = $2;
			}
			my @E_lengths_ar = split(",",$E_lengths);
				
			if (exists $transcripts->{$transcript_name}) {
				$pure_name = $transcript_name;
			} elsif ( (exists $transcripts->{$pure_name}) && ($transcripts_seen{$pure_name} eq $transcript_name) ) {
			
			}  else {
				next;
			}
			
			my $seq = <IN>;
			chomp($seq);
			
			my %ends_seqs;
			
			my $running_total=0;
			for (my $i=0; $i < @E_lengths_ar; $i++) {
				my $exon_seq = substr($seq,$running_total,$E_lengths_ar[$i]);
				$running_total += $E_lengths_ar[$i];
				my $exon_seq_length = (length($exon_seq));
				next if ($exon_seq_length < $MBSL);
				
				if ( $exon_seq_length < $each_side_length ) {
					$ends_seqs{$i}{"B"} = $exon_seq;
					$ends_seqs{$i}{"E"} = $exon_seq;
				} else {
					my $beg_exon_seq = substr($exon_seq,0,($each_side_length));
					my $end_exon_seq = substr($exon_seq,-($each_side_length));
					$ends_seqs{$i}{"B"} = $beg_exon_seq;
					$ends_seqs{$i}{"E"} = $end_exon_seq;
				}
			}
					
			for my $E_Bs (keys %{$transcripts->{$pure_name}}) {
				my ($E, $B);
				if ($E_Bs =~ /E(\d+)B(\d+)/) {
					$E = $1;
					$B = $2;
				}
				
				if (!exists $ends_seqs{$E} || !exists $ends_seqs{$B}) {
					next;
				}
				
				my $bsj_seq = $ends_seqs{$E}{"E"} . $ends_seqs{$B}{"B"};
				
				OUTER:for (my $k = $E; $k < @E_lengths_ar; $k++) {
					if (!exists $ends_seqs{$k}) {
						next;
					}
					for (my $l = $k+1; $l < @E_lengths_ar; $l++) {
						if (!exists $ends_seqs{$l}) {
							next;
						}
						my $fsj_seq = $ends_seqs{$k}{"E"} . $ends_seqs{$l}{"B"};
						if ($bsj_seq eq $fsj_seq) {
							my $fsj = "E$k" . "B$l";
							$bsjs_with_equiv_fsj{$pure_name}{$E_Bs} = $fsj;
							last OUTER;
						}
					}
				}
			}	
		}
	}	
}

sub extract_backsplice_transcripts {
	my ($transcripts, $READ_LENGTHS) = @_;
	open OUT, ">$target_potential_backsplice_fasta" or die "can't open $target_potential_backsplice_fasta\n";
	open OUT_2, ">$target_potential_backsplice_fasta_no_copied_seq" or die "can't open $target_potential_backsplice_fasta_no_copied_seq\n";
	open IN, "<$fa" or die "can't open $fa\n";
	while (<IN>) {
		chomp;
		if ($_ =~ /^>(.*)/) {
			my $transcript_name = $1;
			my ($E_lengths, $pure_name);
			if ($transcript_name =~ m/^.*:\d+-\d+_(.*)_Exon_Lengths_(.*)_Offsets_.*_/) {
					$pure_name = $1;
					$E_lengths = $2;
			}
			my @E_lengths_ar = split(",",$E_lengths);
				
			if (exists $transcripts->{$transcript_name}) {
				$pure_name = $transcript_name;
			} elsif ( (exists $transcripts->{$pure_name}) && ($transcripts_seen{$pure_name} eq $transcript_name) ) {
			
			}  else {
				next;
			}

			my $seq = <IN>;
			chomp($seq);
							
			for my $E_Bs (keys %{$transcripts->{$pure_name}}) {
				my ($e1, $e2);
				if ($E_Bs =~ /E(\d+)B(\d+)/) {
					$e1 = $1;
					$e2 = $2;
				}
				my $running_total = 0;
				my $transcript_sequence;
				if ($e2 != 0) {
					for (my $i=0; $i < $e2; $i++) {
						$running_total += $E_lengths_ar[$i];
					}
				}
				for (my $i = $e2; $i <= $e1; $i++) {
					$transcript_sequence .= substr($seq,$running_total,$E_lengths_ar[$i]);
					$running_total += $E_lengths_ar[$i];
				}
				$running_total -= $E_lengths_ar[$e1];
				
				my $end_anneal_to_beg_sequence;
					
				if ($E_lengths_ar[$e1] < $READ_LENGTHS) {
					$end_anneal_to_beg_sequence = substr($seq,$running_total,$E_lengths_ar[$e1]);
				} else {
					my $end_exon_sequence = substr($seq,$running_total,$E_lengths_ar[$e1]);
					$end_anneal_to_beg_sequence = substr($end_exon_sequence,-($READ_LENGTHS-1));
				}
				
				my $transcript_name = $pure_name . "_" . $E_Bs;
				print OUT_2 ">$transcript_name\n$transcript_sequence\n";
				$transcript_sequence = $end_anneal_to_beg_sequence . $transcript_sequence;
				print OUT ">$transcript_name\n$transcript_sequence\n";				

			}
		}
	}
	close IN;
	close OUT;
	close OUT_2;
}


sub extract_reads {
	my ($outfastq1, $fastq1, $outfastq2, $fastq2, $read_ids) = @_;
	my $reads_list = $mapping_dir . "/read_ids.lst";
	open OUT, ">$reads_list" or die "can't open $reads_list\n";
	for my $read_id (keys %{$read_ids}) {
		print OUT "@" . "$read_id \n";
	}
	close OUT;
	my $command1;
	if ($fastq1 =~ /\.gz$/) {
		$command1 = "gunzip -c $fastq1 | fgrep -A3 --no-group-separator -f $reads_list | gzip > $outfastq1";
	} else {
		$command1 = "fgrep -A3 --no-group-separator -f $reads_list $fastq1 | gzip > $outfastq1";
	}
	
	my $command2;
	if ($fastq2 =~ /\.gz$/) {
		$command2 = "gunzip -c $fastq2 | fgrep -A3 --no-group-separator -f $reads_list | gzip > $outfastq2";
	} else {
		$command2 = "fgrep -A3 --no-group-separator -f $reads_list $fastq2 | gzip > $outfastq2";
	}
	
	my $command = $command1 . "&" . $command2;
	
	system("$command");
}

sub generate_alt_fsj_targets {
	my $each_side_length = $kallisto_k - $MBSL;
	open OUT, ">$fsj_targets_fa" or die "can't open $fsj_targets_fa\n";
	open IN, "<$fa" or die "can't open $fa\n";
	while (<IN>) {
		chomp;
		if ($_ =~ m/^.*:\d+-\d+_(.*)_Exon_Lengths_(.*)_Offsets_.*_/) {
			my $seq = <IN>;
			chomp($seq);
			
			my $pure_name = $1;
			my $E_lengths = $2;
			my @E_lengths_ar = split(",",$E_lengths);
			
			my ($transcript_name) = $_ =~ m/^>(.*)/;
			
			my $raw_name;
			if (exists $transcripts_seen{$transcript_name}) {
				$raw_name = $transcript_name;
			} else {
				$raw_name = $pure_name;
			}
			
			my %ends_seqs;
			
			my $running_total=0;
			for (my $i=0; $i < @E_lengths_ar; $i++) {
				my $exon_seq = substr($seq,$running_total,$E_lengths_ar[$i]);
				$running_total += $E_lengths_ar[$i];
				my $exon_seq_length = (length($exon_seq));
				next if ($exon_seq_length < $MBSL);
				
				if ( $exon_seq_length < $each_side_length ) {
					$ends_seqs{$i}{"B"} = $exon_seq;
					$ends_seqs{$i}{"E"} = $exon_seq;
				} else {
					my $beg_exon_seq = substr($exon_seq,0,($each_side_length));
					my $end_exon_seq = substr($exon_seq,-($each_side_length));
					$ends_seqs{$i}{"B"} = $beg_exon_seq;
					$ends_seqs{$i}{"E"} = $end_exon_seq;
				}
			}
			
			my %outputted_already;
			
			for (my $i=0; $i < @E_lengths_ar; $i++) {
				if (!exists $ends_seqs{$i}) {
					next;
				}
				for (my $j=$i+1; $j < @E_lengths_ar; $j++) {
					if (!exists $ends_seqs{$j}) {
						next;
					}
					
					if ( ($j - $i) > 1) {
						#make sure alternative fsjs have at least one non-skipped exon
						my $at_least_one_non_skipped;
						for (my $k=$i+1; $k < $j; $k++) {
							if (exists $ends_seqs{$k}) {
								$at_least_one_non_skipped = 1;
								last;
							}
						}
						if (!$at_least_one_non_skipped) {
							next;
						}
					} else { #non-alternative
						next;
					}
							
					my $out_header = $raw_name . "_E$i" . "B$j";
					my $out_seq = $ends_seqs{$i}{"E"} . $ends_seqs{$j}{"B"};
					if (exists $outputted_already{$out_header}) {
						next;
					}
					print OUT ">$out_header\n$out_seq\n";
					$outputted_already{$out_header}++;
					
				}
			}
			
			#die;
			
		}
	}
	
	close IN;
	close OUT;

}

sub generate_targets_and_indices {
	
	my $each_side_length = $kallisto_k - $MBSL;
	
	open OUT, ">$all_possible_bsj_targets_fa" or die "can't open $all_possible_bsj_targets_fa\n";
	open OUT2, ">$FirstandLastExons_entities_fa" or die "can't open $FirstandLastExons_entities_fa\n";
	
	open IN, "<$fa" or die "can't open $fa\n";
	while (<IN>) {
		chomp;
		if ($_ =~ m/^.*:\d+-\d+_(.*)_Exon_Lengths_(.*)_Offsets_.*_/) {
			my $seq = <IN>;
			chomp($seq);
			
			my $pure_name = $1;
			my $E_lengths = $2;
			my @E_lengths_ar = split(",",$E_lengths);
			
			my ($transcript_name) = $_ =~ m/^>(.*)/;
			
			my $output_name;
			
			if ( (exists $transcripts_seen{$pure_name}) && ($transcripts_seen{$pure_name} ne $transcript_name) ) {
				$transcripts_seen{$transcript_name} = $transcript_name;
				$pure_name = $transcript_name;
			}
			else {
				$transcripts_seen{$pure_name} = $transcript_name;
			}
			$output_name = $pure_name . "_";
			
			my $tx_num_exons = scalar(@E_lengths_ar);
			if ($tx_num_exons == 1) {
				my $out_header = $output_name . "E0B0";
				print OUT2 ">$out_header\n$seq\n";
			} else {
				my $left_exon_seq= substr($seq, 0, $E_lengths_ar[0]);
				my $right_exon_seq= substr($seq, -($E_lengths_ar[-1]));
				
				my $out_seq = $left_exon_seq. $right_exon_seq;
				my $last_exon = $tx_num_exons - 1;
				my $out_header = $output_name . "E$last_exon" . "B0";
				print OUT2 ">$out_header\n$out_seq\n";
			}
			
			my %ends_seqs;
			
			my $running_total=0;
			for (my $i=0; $i < @E_lengths_ar; $i++) {
				my $exon_seq = substr($seq,$running_total,$E_lengths_ar[$i]);
				$running_total += $E_lengths_ar[$i];
				my $exon_seq_length = (length($exon_seq));
				next if ($exon_seq_length < $MBSL);
				
				if ( $exon_seq_length < $each_side_length ) {
					$ends_seqs{$i}{"B"} = $exon_seq;
					$ends_seqs{$i}{"E"} = $exon_seq;
				} else {
					my $beg_exon_seq = substr($exon_seq,0,($each_side_length));
					my $end_exon_seq = substr($exon_seq,-($each_side_length));
					$ends_seqs{$i}{"B"} = $beg_exon_seq;
					$ends_seqs{$i}{"E"} = $end_exon_seq;
				}
			}
			
			for (my $i=0; $i < @E_lengths_ar; $i++) {
				if (!exists $ends_seqs{$i}) {
					next;
				}
				for (my $j=$i; $j < @E_lengths_ar; $j++) {
					if (!exists $ends_seqs{$j}) {
						next;
					}
					my $out_header = $output_name . "E$j" . "B$i";
					my $out_seq = $ends_seqs{$j}{"E"} . $ends_seqs{$i}{"B"};
					print OUT ">$out_header\n$out_seq\n";
					
				}
			}
			
			#die;
			
		}
	}
	close IN;
	close OUT;
	close OUT2;
	
	system("$kallisto index -i $all_possible_bsj_targets_fa_index -k $kallisto_k $all_possible_bsj_targets_fa");
	system("$kallisto index -i $FirstandLastExons_entities_fa_index -k  $final_kallisto_k $FirstandLastExons_entities_fa");
}

sub rev_compl {
	my $DNA_seq = shift;
	$DNA_seq = reverse($DNA_seq);
	$DNA_seq =~ tr/ATCG/TAGC/;
	return $DNA_seq;
}
