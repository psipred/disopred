#!/usr/bin/perl -w
##
## Output formatting for disopred3
##

use strict;
use File::Basename;
use Cwd qw(abs_path);

# format_protein_binding_predictions($diso3_fn, $svc_fn, $pb_fn);
format_protein_binding_predictions($ARGV[0], $ARGV[1], $ARGV[2]);


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
