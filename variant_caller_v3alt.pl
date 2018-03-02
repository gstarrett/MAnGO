#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::DB::HTS;
use Getopt::Long;

#my $infile = shift; # Input fastq files for alignment
my $index = shift;
my $ref = shift; # Reference fasta
my $range = shift;
my $trans= shift;

# structure of output files
# uid	npos	nref	nalt	nmut	codon	trint apos aalt	amut	aachange	count	freq	Escore	file

my (@summary, @prefilter, @out, %readLenHash, %posHash, %readHash, %baseHash, %simpleReads,%uniq);
my ($bamPath,$samPath,$sortPath);

open (INDEX, "< $index");
my @indexFiles = <INDEX>;

for my $infile (@indexFiles) {
	$sortPath = $infile;
	### high level API
	my $sam = Bio::DB::HTS->new(-bam => $sortPath, -fasta=> $ref, -force_refseq => 1) or die "cannot load bam file";
	my @targets = $sam->seq_ids;

	my $finalLen = 0;

	my @r1 = split(",",$range);
	for (@r1) {
	 my @r2 = split("-",$_);
	 my $delLen = 0;
	 my $snp_caller = sub {
	 			 my ($seqid,$pos,$p) = @_;
	 			 my $refbase = uc($sam->segment($seqid,$pos,$pos)->dna);
				 if ($pos >= $r2[0] && $pos <= $r2[1]) {
					 if ($delLen > 0) {
						 $delLen--;
						 next;
					 }
					 $simpleReads{"ref"} .= $refbase;
	 			 for my $pileup (@$p) {
	 					 my $b = $pileup->alignment;
	 					 my $inBool = $pileup->indel;
	 					 my $qreadname = $b->query->name;
	 					 my $qbase = uc(substr($b->qseq,$pileup->qpos,1));
	 					 if ($inBool > 0) {
							 $qbase = "+";
	 						 $readHash{$qreadname}++;
	 					 } elsif ($inBool < 0) {
	 						 $readHash{$qreadname}++;
							 $qbase = "-" x abs($inBool);
							 $delLen = $inBool;
	 					 } else {
	 						 my $qscore = $b->qscore->[$pileup->qpos];
	 						 $qbase = "N" unless $qscore >= 30;
	 				 	}
	 					$simpleReads{$qreadname} .= $qbase;
	 			 }
			 } else {
				 $simpleReads{$qreadname} .= "N";
			 }
	 };
	 $sam->pileup("$targets[0]:$r2[0]-$r2[1]",$snp_caller);
	 	$finalLen += $r2[1]-$r2[0]+1;
	}
}

for my $key (keys %simpleReads) {
	my $simpleRead = $simpleReads{$key};
	next unless length($simpleRead) == $finalLen;
	if (exists $uniq{$simpleRead}) {
		${$uniq{$simpleRead}}[0]++;
	} else {
		$uniq{$simpleRead} = [1,&generateUID;];
	}
}
open(OUT, "> $infile.counts.uniq.txt");
for my $uniqSeq (keys %uniq) {
	my @refArray = split('', $simpleReads{"ref"});
	my @seqArray = split('', $uniqSeq);
	for (my $i = 1; $i < $#seqArray-1; $i++) {
		my $oline;
		my $ref = $refArray[$i];
		my $alt = $seqArray[$i];
		if ($alt ne $ref) {
			my $trint = join("", $seqArray[$i-1], "x", $seqArray[$i+1]);
			my $apos = ;
			my $aref = ;
			my $aalt = ;
			my $amut = "$aref>$aalt";
			my $aachange;
			if ($aref ne $aalt) {
				$aachange = "Y";
			} else {
				$aachange = "N";
			}
			my $count = $uniq{$uniqSeq};
			my $freq = $count/$total;
			my $enrich = $freq/$reffreq;
			$oline = join("\t", $uid, $i, $ref, $alt, "$ref>$alt", $codon, $trint, $apos, $aalt,	$amut,	$aachange, $count, $freq, $enrich, $infile);
		}
	}
}
close(OUT);

open(OUT, "> $infile.counts.txt");

foreach (@out) {
	print OUT $_, "\n";
}

close(OUT);

#############
#SUBROUTINES#
#############


sub timeStamp {
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	my $year = 1900 + $yearOffset;
	my $theTime = $hour . $minute . $second . "_" . $year . "-" . $month . "-" . $dayOfMonth;
	return($theTime)
}

sub translatorMain {
		my $dna = shift;
		$dna =~ s/-//g;
		my @proteins;
		my $protein;
		for(my $i=shift; $i < (length($dna) - 2) ; $i += 3) {
			my $subseq = substr($dna,$i,3);
			$protein .= &translator($subseq);
		}
		push (@proteins, $protein);
		return @proteins;
}

sub translator {
	my %daRules = (
	'TCA' => 'S',    # Serine
	'TCC' => 'S',    # Serine
	'TCG' => 'S',    # Serine
	'TCT' => 'S',    # Serine
	'TTC' => 'F',    # Phenylalanine
	'TTT' => 'F',    # Phenylalanine
	'TTA' => 'L',    # Leucine
	'TTG' => 'L',    # Leucine
	'TAC' => 'Y',    # Tyrosine
	'TAT' => 'Y',    # Tyrosine
	'TAA' => '*',    # Stop
	'TAG' => '*',    # Stop
	'TGC' => 'C',    # Cysteine
	'TGT' => 'C',    # Cysteine
	'TGA' => '*',    # Stop
	'TGG' => 'W',    # Tryptophan
	'CTA' => 'L',    # Leucine
	'CTC' => 'L',    # Leucine
	'CTG' => 'L',    # Leucine
	'CTT' => 'L',    # Leucine
	'CCA' => 'P',    # Proline
	'CCC' => 'P',    # Proline
	'CCG' => 'P',    # Proline
	'CCT' => 'P',    # Proline
	'CAC' => 'H',    # Histidine
	'CAT' => 'H',    # Histidine
	'CAA' => 'Q',    # Glutamine
	'CAG' => 'Q',    # Glutamine
	'CGA' => 'R',    # Arginine
	'CGC' => 'R',    # Arginine
	'CGG' => 'R',    # Arginine
	'CGT' => 'R',    # Arginine
	'ATA' => 'I',    # Isoleucine
	'ATC' => 'I',    # Isoleucine
	'ATT' => 'I',    # Isoleucine
	'ATG' => 'M',    # Methionine
	'ACA' => 'T',    # Threonine
	'ACC' => 'T',    # Threonine
	'ACG' => 'T',    # Threonine
	'ACT' => 'T',    # Threonine
	'AAC' => 'N',    # Asparagine
	'AAT' => 'N',    # Asparagine
	'AAA' => 'K',    # Lysine
	'AAG' => 'K',    # Lysine
	'AGC' => 'S',    # Serine
	'AGT' => 'S',    # Serine
	'AGA' => 'R',    # Arginine
	'AGG' => 'R',    # Arginine
	'GTA' => 'V',    # Valine
	'GTC' => 'V',    # Valine
	'GTG' => 'V',    # Valine
	'GTT' => 'V',    # Valine
	'GCA' => 'A',    # Alanine
	'GCC' => 'A',    # Alanine
	'GCG' => 'A',    # Alanine
	'GCT' => 'A',    # Alanine
	'GAC' => 'D',    # Aspartic Acid
	'GAT' => 'D',    # Aspartic Acid
	'GAA' => 'E',    # Glutamic Acid
	'GAG' => 'E',    # Glutamic Acid
	'GGA' => 'G',    # Glycine
	'GGC' => 'G',    # Glycine
	'GGG' => 'G',    # Glycine
	'GGT' => 'G',    # Glycine
	);

	my $codon = shift;
	if ($codon =~ /N/) {
		return "X";
	} elsif ($codon =~ /[+-]/) {
		return "#";
	} elsif (exists $daRules{$codon}) {
		return $daRules{$codon};
	} else {die "Bad move Anton: $codon\n";}
}
sub translation_align {
	my $seq1 = shift;
	my $seq2 = shift;
	my $seqLength = length($seq1);
	#print $seqLength;
	my $finalseq;

	my $i=0;
	until ($i == $seqLength){
		my $subseq1 = substr($seq1,$i,1);
		my $subseq2 = substr($seq2,$i,1);
		if ($subseq2 =~ /$subseq1/) {
			$finalseq .= '.';
		} else {
			$finalseq .= "$subseq2";
		}
		$i++;
	}
}
sub generateUID {
	my $uid;
	my @hexLib = (0..9,"a".."z");
	for (my $var = 0; $var < 12; $var++) {
		$uid .= int(rand($#hexLib));
	}
	return($uid);
}
