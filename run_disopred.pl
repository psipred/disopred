#!/usr/bin/perl -w
##
## *********************************************************************
## *    DISOPRED3 - PROTEIN INTRINSIC DISORDER PREDICTION PROGRAM BY   *
## *                  D.T. JONES, D. COZZETTO & J.J. WARD              *
## *********************************************************************
##
## Please refer to the license agreement for permission of use and re-distribution
##

use strict;
use File::Basename;
use Cwd qw(abs_path);


## IMPORTANT: Set the paths to folder with the NCBI executables and to the
## sequence database
my $NCBI_DIR = "/scratch0/NOT_BACKED_UP/dbuchan/Applications/blast-2.2.26/bin/";
my $SEQ_DB = "/scratch0/NOT_BACKED_UP/dbuchan/uniref/uniref_test_db/uniref_test.fasta";

## IMPORTANT: Changing these flags will alter the behaviour of blastpgp
## You may want to use -a n to speed-up the search using n processors, if available
my $PSIBLAST_PAR = "-a 1 -b 0 -j 3 -h 0.001";

## IMPORTANT: Moving the bin/, data, or dso_lib directories to a different location will cause the programs
## to crash, unless you change the variables below accordingly
my $dir = dirname $0;
my $EXE_DIR = abs_path(join '/', $dir, "bin"); # the path of the bin directory
my $DATA_DIR = abs_path(join '/', $dir,"data"); # the path of the data directory
$ENV{DSO_LIB_PATH} = join '/', abs_path($dir), "dso_lib/"; # the path of the library directory used by the nearest neighbour classifier
exists $ENV{DSO_LIB_PATH} or die "[$0] ERROR: DSO_LIB_PATH environmental variable not set\n";

my $DISO2_FPR = 5; # the adjustable DISOPRED2 false positive rate, represented as an integer between 1 and 10

# check whether the command line contains only one argument
die "[$0] ERROR: The command line should only contain the path to the input sequence file\n" if ( scalar @ARGV != 1);

# die if input file does not exist or is not a text file
die "[$0] ERROR: Input file $ARGV[0] does not exist\n"  if !-e $ARGV[0];
die "[$0] ERROR: Input file $ARGV[0] does not look like a text file\n"  if !-T $ARGV[0];

my $fasta_fn = abs_path($ARGV[0]);
my ($out_dir, $base) = (dirname($fasta_fn), basename($fasta_fn));

$base =~ s/\.fa(sta)?$//;
my $host_id = `hostid`;
chomp $host_id;
my $tmp_base = join '_', $base , $$, $host_id;

my ($hits_file, $chk_file) = map {my $name = join '.', $tmp_base, $_; join '/', $out_dir, $name } ("blast", "chk");

print "Running PSI-BLAST search ...\n\n";
# run psiblast
my $args = join ' ', $NCBI_DIR."blastpgp", "-i", $fasta_fn, "-d", $SEQ_DB, $PSIBLAST_PAR, "-C", $chk_file, "-o", $hits_file, "\n";
system($args) == 0 or die "[$0] ERROR: $args failed: $?\n";

print "Generating PSSM ...\n\n";

$args = join ' ', "echo", $chk_file, ">", $tmp_base.".pn", "\n";
system($args) == 0 or die "[$0] ERROR: $args failed: $?\n";

$args = join ' ', "echo", $fasta_fn, ">", $tmp_base.".sn", "\n";
system($args) == 0 or die "[$0] ERROR: $args failed: $?\n";

$args = join ' ', $NCBI_DIR."makemat", "-P", $tmp_base, "\n";
system($args) == 0 or die "[$0] ERROR: $args failed: $?\n";

my $mtx_fn = join '/', $out_dir, $tmp_base.".mtx";
die "[$0] ERROR: Couldn't find the mtx file $mtx_fn\n" if !-e $mtx_fn;

my @exts = ("diso", "diso2", "nndiso", "dnb", "diso", "in_svm_dat", "out_svm_dat", "pbdat");
my ($diso_fn, $diso2_fn, $nndiso_fn, $dnb_fn, $diso3_fn, $dat_fn, $svc_fn, $pb_fn) = map { abs_path(join '/', $out_dir, $base.".$_") } @exts;

print "Predicting disorder with DISOPRED2 ...\n\n";
$args = join ' ', "$EXE_DIR/disopred2", join('/', $out_dir, $base), $mtx_fn, "$DATA_DIR/", $DISO2_FPR, "\n";
system($args) == 0 or die "[$0] ERROR: $args failed. Please report error to psipred\@cs.ucl.ac.uk\n";

$args = join ' ', "mv", $diso_fn, $diso2_fn, "\n";
system($args) == 0 or die "[$0] ERROR: $args failed. Please report error to psipred\@cs.ucl.ac.uk";

print "Running neural network classifier ...\n\n";
$args = join ' ', "$EXE_DIR/diso_neu_net", "$DATA_DIR/weights.dat.nmr_nonpdb", $mtx_fn, ">", $nndiso_fn, "\n";
system($args) == 0 or die "[$0] ERROR: $args failed. Please report error to psipred\@cs.ucl.ac.uk";

print "Running nearest neighbour classifier ...\n\n";
$args = join ' ', "$EXE_DIR/diso_neighb", $mtx_fn, "$DATA_DIR/dso.lst", ">", $dnb_fn, "\n";
system($args) == 0 or die "[$0] ERROR: $args failed. Please report error to psipred\@cs.ucl.ac.uk\n";

print "Combining disordered residue predictions ...\n\n";
$args = join ' ', "$EXE_DIR/combine", "$DATA_DIR/weights_comb.dat", $diso2_fn, $nndiso_fn, $dnb_fn, ">", $diso3_fn, "\n";
system($args) == 0 or die "[$0] ERROR: $args failed. Please report error to psipred\@cs.ucl.ac.uk";

my ($seq, $idr_data) = parse_disopred3_file($diso3_fn);
# parse mtx file and extract the profile data
my $profile = get_lines_from_mtx_file($mtx_fn);
my $feat_vecs = make_vectors($profile, $seq, $idr_data, 15);

open(DAT, '>', $dat_fn) or die "[$0] ERROR: Couldn't open output file $dat_fn\n";
print DAT join "\n", @{$feat_vecs}, '';
close DAT;

print "Predicting protein binding residues within disordered regions ...\n\n";
$args = join ' ', "$EXE_DIR/svm-predict", "-b 1 -q", $dat_fn, "$DATA_DIR/ProtBind_IDR_Model.dat", $svc_fn;
system($args) == 0 or die "[$0] ERROR: $args failed. Please report error to psipred\@cs.ucl.ac.uk\n"; #system $args failed: $?\n";

# merge residue level predictions of disorder and protein binding
format_protein_binding_predictions($diso3_fn, $svc_fn, $pb_fn);

# Remove temporary files
print "Cleaning up ...\n\n";
$args = join ' ', "rm -f", $hits_file, $chk_file, "error.log", $mtx_fn, <$tmp_base*>, glob("$out_dir/*horiz_d"), $diso2_fn, $nndiso_fn, $dnb_fn, $dat_fn, $svc_fn, "\n";
system($args) == 0 or die "[$0] ERROR: $args failed: $?\n";

print join "\n\n", "Finished", "Disordered residue predictions in $diso3_fn", "Protein binding disordered residue predictions in $pb_fn", '';
0;

# Parse disordered residue predictions and obtain positional information about intrinsically disordered regions
sub parse_disopred3_file {
	my $pred_fn = $_[0];
	open(DISO, $pred_fn) or die "[$0] ERROR: Couldn't open Disopred output file $pred_fn\n";

	my (@idr_id, @aa, %starts, %ends, %lengths) = ((), (),(), (), ());
	my ($cur_seg, $cur_length, $last_pos) = (0, 0, -1);
	while (<DISO>) {
		if ($_ =~ m/^\s*(\d+)\s([A-Z])\s[\.\*]\s(\S+)/ ) {
			push @aa, $2;

			if ($3 >= 0.5) {
				$cur_seg++ if $1 != $last_pos+1;
				$cur_length++;
				push @idr_id, $cur_seg;
				$starts{$cur_seg} = $1 if $1 != $last_pos+1;
				if (eof DISO) {
					$ends{$cur_seg} = $1;
					$lengths{$cur_seg} = $cur_length;
					$cur_length = 0
				}
				$last_pos = $1;
			}

			else {
				push @idr_id, 0;
				if ($1 == $last_pos+1) {
					$ends{$cur_seg} = $last_pos;
					$lengths{$cur_seg} = $cur_length;
					$cur_length = 0
				}
			}
		}
	}

	close DISO;
	my @data = map { $_  ? [sprintf("%.6f", log(1 + $lengths{$_})), sprintf("%.6f", $starts{$_}/(scalar @aa)), sprintf("%.6f", $ends{$_}/(scalar @aa))] : [0,0,0] } @idr_id;
	scalar @aa == scalar @data or die "[$0] ERROR: Different number of amino acids and disorder region data vectors from $pred_fn\n";
	return (\@aa, \@data)
}


sub get_lines_from_mtx_file {
	my $mtx_file = $_[0];

	open(MTX, $mtx_file) or die "[$0] ERROR: Couldn't open $mtx_file makemat output file\n";
	my @lines = <MTX>;
	close MTX;
	chop @lines;
	scalar @lines == $lines[0] + 14 or die "[$0] ERROR: Unexpected number of lines in $mtx_file\n";

	my $par = get_linear_scaling_params();
	my @profile = ();
	foreach my $pos (1 .. $lines[0]) {
		print join ' ', "Undefined line", $pos, $pos + 13, "\n" if !defined $lines[ $pos+13 ];
		my @data = split /\s+/, $lines[$pos+13];
		# extract current profile data for standard amino acids
		my @pos = (1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22);
		my $par = get_linear_scaling_params();

		# linearly scale PSSM values based on the range of scores observed while training
		my @scaled_data = map { $_ ne 21 ?
			sprintf "%.6f", ($data[$_] - $$par{$_}{'min'})/($$par{$_}{'max'} - $$par{$_}{'min'} ) :
			sprintf "%.6f", $data[$_] - $$par{$_}{'min'} } @pos;
		push @profile,  [ @scaled_data ] ;
	}

	return \@profile
}

# Read in the maximum and minimum PSSM scores observed in the training data after 3 iterations of PSIBLAST
sub get_linear_scaling_params {
	my %param_linear_scaling = ();

        $param_linear_scaling{1}{'min'}  =  -956;
        $param_linear_scaling{1}{'max'}  =   734;
        $param_linear_scaling{3}{'min'}  = -1021;
        $param_linear_scaling{3}{'max'}  =  1353;
        $param_linear_scaling{4}{'min'}  = -1120;
        $param_linear_scaling{4}{'max'}  =   899;
        $param_linear_scaling{5}{'min'}  = -1063;
        $param_linear_scaling{5}{'max'}  =   827;
        $param_linear_scaling{6}{'min'}  = -1002;
        $param_linear_scaling{6}{'max'}  =  1071;
        $param_linear_scaling{7}{'min'}  = -1005;
        $param_linear_scaling{7}{'max'}  =   808;
        $param_linear_scaling{8}{'min'}  =  -995;
        $param_linear_scaling{8}{'max'}  =  1287;
        $param_linear_scaling{9}{'min'}  = -1047;
        $param_linear_scaling{9}{'max'}  =   866;
        $param_linear_scaling{10}{'min'} = -1001;
        $param_linear_scaling{10}{'max'} =   858;
        $param_linear_scaling{11}{'min'} = -1006;
        $param_linear_scaling{11}{'max'} =   719;
        $param_linear_scaling{12}{'min'} =  -954;
        $param_linear_scaling{12}{'max'} =  1204;
        $param_linear_scaling{13}{'min'} = -1070;
        $param_linear_scaling{13}{'max'} =   922;
        $param_linear_scaling{14}{'min'} = -1081;
        $param_linear_scaling{14}{'max'} =   909;
        $param_linear_scaling{15}{'min'} =  -986;
        $param_linear_scaling{15}{'max'} =   967;
        $param_linear_scaling{16}{'min'} = -1039;
        $param_linear_scaling{16}{'max'} =   936;
        $param_linear_scaling{17}{'min'} =  -975;
        $param_linear_scaling{17}{'max'} =   779;
        $param_linear_scaling{18}{'min'} =  -945;
        $param_linear_scaling{18}{'max'} =   827;
        $param_linear_scaling{19}{'min'} =  -998;
        $param_linear_scaling{19}{'max'} =   780;
        $param_linear_scaling{21}{'min'} =  -100;
        $param_linear_scaling{21}{'max'} =  -100;
        $param_linear_scaling{22}{'min'} =  -995;
        $param_linear_scaling{22}{'max'} =  1107;
	return \%param_linear_scaling
}

sub make_vectors {
	my ($prf, $res, $length_pos_data, $win_size) = @_;

	my $n_col = scalar @$prf;
	$n_col == scalar @$length_pos_data or die "[$0] ERROR: Different numbers of elements in the profile data structure and the array of disordered region lengths\n";
	my @lines = ();

	for (my $i = 0; $i < $n_col; $i++) {
		my $flag = 0;
		my ($start, $end) = ($i - ($win_size-1)/2, $i + ($win_size-1)/2);
		my ($first_label, $last_label) = (1, 20*$win_size);

		while ($start < 0) {
			$flag = 1 if !$flag;
			$start++;
			$first_label += 20;
		}

		while ($end >= $n_col ) {
			$flag = 1 if !$flag;
			$end--;
			$last_label -= 20;
		}

		my @f_values = ();
		foreach my $el ($start..$end) {
			push @f_values, @{$$prf[$el]}
		}

		my @f_indexes = $first_label..$last_label;

		# append to the scaled profile data the flag for windows exceeding the input sequence, the positional information of
		# any predicteddisordered region and the amino acid composition in the current window
		push @f_values, $flag, @{$$length_pos_data[$i]};
		push @f_indexes, 20*$win_size+1 .. 20*$win_size+4;

		my @alphabet = ("A", "C", "D", "E", "F", "G", "H", "K", "I", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y");
		my $cur_seq = join '', @$res[$start..$end];
		my $l = length $cur_seq;
		my @aa_comp = map { sprintf "%.6f", ($cur_seq =~ s/$_/$_/g)/$l } @alphabet;
		push @f_values, @aa_comp;
		push @f_indexes,  20*$win_size+5 .. 20*$win_size+24;

		@f_indexes == @f_values or die "[$0] ERROR: Different number of feature values (", scalar @f_values , ") and labels (", scalar @f_indexes ,")\n";

		my @data = map { $f_values[$_] > 0 ? ( join ':', $f_indexes[$_], $f_values[$_] ) : () } 0..scalar(@f_indexes)-1;
		my $sub_seq = join '', @$res[$start..$end];
		my $cur_line = join " ", 0, @data;
		$cur_line = join " # ", $cur_line, $sub_seq;
		push @lines, $cur_line
	}
	return \@lines
}



sub format_protein_binding_predictions {
	my ($idr_fn, $pb_fn, $out_fn) = @_;
	open(DISO, $idr_fn) or die "[$0] ERROR: Couldn't open Disopred output file $idr_fn\n\n";

	my (@idr_scores, @aa) = ((),());
	while (<DISO>) {
		if ($_ =~ m/^\s*\d+\s([A-Z])\s[\.\*]\s(\S+)/ ) {
			push @aa, $1;
			push @idr_scores, $2
		}
	}
	close DISO;
	scalar @aa == scalar @idr_scores or die "[$0] ERROR: Uneven number of amino acids and disorder confidence scores from $idr_fn\n\n";

	open(PB, $pb_fn) or die "[$0] ERROR: Couldn't open Disopred output file $pb_fn\n";

	my (@conf_scores, @pred_classes) = ((), ());
	while (<PB>) {
		chop;
		my @tokens = split /\s+/;
		scalar @tokens == 3 or die "[$0] ERROR: Unexpected number of fields (not 3) at line\n$_\nin $pb_fn\n\n";
		my $k = 1;
		if ($tokens[0] eq "labels") {
			shift @tokens;
			scalar @tokens == 2 or die "[$0] ERROR: Unexpected number of labels in $pb_fn\n\n";
			$k = 2 if $tokens[0] != 1;
		}

		else  {
			push @conf_scores, $tokens[$k]
		}
	}
	close PB;

	scalar @conf_scores == scalar @aa or die "[$0] ERROR: Uneven number of class assignments and amino acids from $pb_fn and $idr_fn\n\n";

	my @lines = ();
	for (my $i = 0; $i < scalar @aa; $i++) {
		if ($idr_scores[$i] >= 0.5) {
			my $cur_state = $conf_scores[$i] >= 0.5 ? "^" : "-";
			push @lines, sprintf("%5d %s %s %4.2f", $i+1, $aa[$i], $cur_state, $conf_scores[$i])
		}
		else {
			push @lines, sprintf("%5d %s %s %4s", $i+1, $aa[$i], ".", "NA")
		}
	}

	open(OUT, '>', $out_fn) or die "[$0] ERROR: Couldn't open output file $out_fn\n";
	print OUT "#                  ----- DISOPRED version 3.1 -----\n";
	print OUT "#     Protein binding site prediction within disordered regions\n";
	print OUT "#   Protein-binding disordered residues are marked with carets (^)\n";
	print OUT "# Disordered residues not binding proteins are marked with dashes (-)\n";
	print OUT "#            Ordered amino acids are marked with dots (.)\n";
	print OUT join "\n", @lines, '';
	close OUT
}
