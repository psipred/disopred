#!/usr/bin/perl -w

#
# Small script to abstract some components for the website
#

use strict;
use File::Basename;
use Cwd qw(abs_path);

my $diso3_fn = $ARGV[0];
my $mtx_fn = $ARGV[1];
my $dat_fn = $ARGV[2];

my ($seq, $idr_data) = parse_disopred3_file($diso3_fn);
# parse mtx file and extract the profile data
my $profile = get_lines_from_mtx_file($mtx_fn);
my $feat_vecs = make_vectors($profile, $seq, $idr_data, 15);

open(DAT, '>', $dat_fn) or die "[$0] ERROR: Couldn't open output file $dat_fn\n";
print DAT join "\n", @{$feat_vecs}, '';
close DAT;

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
