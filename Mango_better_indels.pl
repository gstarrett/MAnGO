#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::DB::HTS;
use Getopt::Long;

#my $infile = shift; # Input fastq files for alignment
my $index = shift;
my $ref = shift; # Reference fasta
my $range = shift;
my $prefix = shift;
my $codonStart = shift;
my %lengthHash;

my $cdhit = "/Applications/bioinf_tools/cdhit/cd-hit-est";
my $blat = "/Applications/bioinf_tools/blat";

# structure of output files
# uid	npos	nref	nalt	nmut	codon	trint apos aalt	amut	aachange	count	freq	Escore	file

my (%uniq,%uids,$insertSize,%fileHash);
my ($bamPath,$samPath,$sortPath);
my $min = 999999999999;
my $max = 0;
#my $finalLen = 0;
my %simpleReads;

open (INDEX, "< $index");
my @indexFiles = <INDEX>;
close(INDEX);
my $insertStart = 0;
my $insertEnd = 0;
my @r1 = split(",",$range);
for (@r1) {
 my @r2 = split("-",$_);
 $insertStart = $r2[1] if $insertStart == 0;
 $min = $r2[0] if $r2[0] < $min;
 $max = $r2[1] if $r2[1] > $max;
 $insertEnd = $r2[0] if $r2[0] > $insertEnd;
}
open (FASTA, "> $prefix.fa");
for (my $n = 0; $n <= $#indexFiles; $n++) {
	%simpleReads = ();
	my $infile = $indexFiles[$n];
	chomp($infile);
	print "Processing $infile...\n";
	$sortPath = $infile;
	### high level API
	my $sam = Bio::DB::HTS->new(-bam => $sortPath, -fasta=> $ref, -force_refseq => 1, -autoindex => 1) or die "cannot load bam file";
	my @targets = $sam->seq_ids;
	 my $delLen = 0;
	 my $snp_caller = sub {
	 			 my ($seqid,$pos,$p) = @_;
	 			 my $refbase = uc($sam->segment($seqid,$pos,$pos)->dna);
				 if (($pos >= $min && $pos <= $insertStart) || ($pos >= $insertEnd && $pos <= $max)) {
					 $simpleReads{"ref"} .= $refbase;
	 			 	 for my $pileup (@$p) {
						 my $b = $pileup->alignment;
						 my $qreadname = $b->query->name;
						 if ($delLen > 0) {
							 $delLen--;
							 next;
						 }
	 					 my $inBool = $pileup->indel;
	 					 my $qbase = uc(substr($b->qseq,$pileup->qpos,$inBool));
	 					 # if ($inBool > 0) {
							#  $qbase .= "+";
	 					 # } elsif
             if ($inBool < 0) {
							 $qbase = "-" x abs($inBool);
							 $delLen = abs($inBool);
	 					 } else {
	 						 my $qscore = $b->qscore->[$pileup->qpos];
	 						 if ($qscore <25 && $qbase ne $refbase) {
                 print "Found low quality ambiguous base: $pos:$refbase>$qbase:$qscore\n";
                 $qbase = "N";
               }
	 				 	}
	 					$simpleReads{$qreadname} .= $qbase;
	 			 }
			 } elsif ($pos > $insertStart && $pos < $insertEnd) {
         print "Inputting $pos:N for insert\n";
				 $simpleReads{"ref"} .= "N";
				 for my $pileup (@$p) {
					 my $b = $pileup->alignment;
					 my $qreadname = $b->query->name;
				 	 $simpleReads{$qreadname} .= "N";
			 	 }
			 }
	 };
	 print "\tFinding reads overlapping targets...\n";
	 $sam->pileup("$targets[0]:$min-$max",$snp_caller);

	for my $key (keys %simpleReads) {
    next if $key eq "ref";
		my $simpleRead = $simpleReads{$key};
		#$finalLen = $max - $min;
		#print length($simpleRead) . "\t" . $finalLen . "\t" . length($simpleReads{"ref"}) . "\n";
    if (exists $lengthHash{length($simpleRead)}) {
      $lengthHash{length($simpleRead)} = 1;
    } else {
      $lengthHash{length($simpleRead)}++;
    }
		next unless length($simpleRead) == length($simpleReads{"ref"});
		if (exists $uniq{$simpleRead}) {
			my $existUID = $uniq{$simpleRead};
			if (exists ${$uids{$existUID}}{$infile}) {
        my $current = ${$uids{$existUID}}{$infile};
				${$uids{$existUID}}{$infile}++;
        print "Found existing UID in $infile:$existUID:$current\n";
			} else {
				${$uids{$existUID}}{$infile} = 1;
        print "Found first presence of UID:$existUID in $infile\n";
			}
		} else {
      my $newUID = &generateUID($n);
      print "Found new sequence in $infile... assigning UID:$newUID\n";
			$uniq{$simpleRead} = $newUID;
			my %internal = ($infile => 1);
			$uids{$newUID} = \%internal;
      print FASTA ">" . $newUID . "\n$simpleRead\n"
		}
		if (exists $fileHash{$infile}) {
			$fileHash{$infile}++;
		} else {
			$fileHash{$infile} = 1;
		}
	}
}
close(FASTA);
print "Completed finding reads!\n";
open(OUT, "> $prefix.counts.uniq.txt");
my $seqCount = scalar keys %uniq;
print "Calling variants within $seqCount unique sequences...\n";
for my $uniqSeq (keys %uniq) {
	my $uid = $uniq{$uniqSeq};
	for my $infile (keys %{$uids{$uid}}) {
		my $total = $fileHash{$infile};
		my @refArray = split('', $simpleReads{"ref"});
		my @seqArray = split('', $uniqSeq);
		for (my $i = 1; $i < $#seqArray-1; $i++) {
			my ($pos1, $pos2, $pos3, $refpos1, $refpos2, $refpos3, $refcodon, $codon, $aref, $aalt);
			my $oline;
			my $ref = $refArray[$i];
			my $alt = $seqArray[$i];
			next if $alt eq "N";
			my $codonPos = ($i+$codonStart-1)%3;
			my $apos = int(($i+$codonStart-1)/3);

			if ($alt ne $ref) {
				if ($codonPos == 0) {
					$pos1 = $alt;
					$pos2 = lc($seqArray[$i+1]);
					$pos3 = lc($seqArray[$i+2]);
					$refpos1 = $ref;
					$refpos2 = lc($refArray[$i+1]);
					$refpos3 = lc($refArray[$i+2]);
					$codon = join("", $pos1, $pos2, $pos3);
					$refcodon = join("", $refpos1, $refpos2, $refpos3);
					$aref = &translator(uc($refcodon));
					$aalt = &translator(uc($codon));
				} elsif ($codonPos == 1) {
					$pos1 = lc($seqArray[$i-1]);
					$pos2 = $alt;
					$pos3 = lc($seqArray[$i+1]);
					$refpos1 = lc($refArray[$i-1]);
					$refpos2 = $ref;
					$refpos3 = lc($refArray[$i+1]);
					$codon = join("", $pos1, $pos2, $pos3);
					$refcodon = join("", $refpos1, $refpos2, $refpos3);
					$aref = &translator(uc($refcodon));
					$aalt = &translator(uc($codon));
				} elsif ($codonPos == 2) {
					$pos1 = lc($seqArray[$i-2]);
					$pos2 = lc($seqArray[$i-1]);
					$pos3 = $alt;
					$refpos1 = lc($refArray[$i-2]);
					$refpos2 = lc($refArray[$i-1]);
					$refpos3 = $ref;
					$codon = join("", $pos1, $pos2, $pos3);
					$refcodon = join("", $refpos1, $refpos2, $refpos3);
					$aref = &translator(uc($refcodon));
					$aalt = &translator(uc($codon));
				}
				my $trint = join("", $seqArray[$i-1], "x", $seqArray[$i+1]);
				my $amut = "$aref>$aalt";
				my $aachange;
				if ($aref ne $aalt) {
					$aachange = "Y";
				} else {
					$aachange = "N";
				}
				my $count = ${$uids{$uid}}{$infile};
				my $freq = $count/$total;
				#my $enrich = $freq/$reffreq;
				$oline = join("\t", $uid, $i+$min, $ref, $alt, "$ref>$alt", $codon, $trint, $apos, $aref, $aalt,	$amut, $aachange, $count, $freq, $infile);
				print OUT $oline,"\n";
			}
		}
	}
}
close(OUT);
print "Analysis complete!\n";
print Dumper(%lengthHash);
#############
#SUBROUTINES#
#############


sub timeStamp {
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	my $year = 1900 + $yearOffset;
	my $theTime = $hour . $minute . $second . "_" . $year . "-" . $month . "-" . $dayOfMonth;
	return($theTime)
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
	if ($codon =~ /N/ || length($codon) != 3) {
		return "X";
	} elsif ($codon =~ /[+-]/) {
		return "^";
	} elsif (exists $daRules{$codon}) {
		return $daRules{$codon};
	} else {die "Bad move Anton: $codon\n";}
}

sub generateUID {
	my $uid;
	my @hexLib = (0..9,"a".."z");
	for (my $var = 0; $var < 12; $var++) {
		$uid .= $hexLib[int(rand($#hexLib))];
	}
	return($uid);
}

# sub blatResults {
#   # goal identify rearrangments and report in a logical fashion
#   # bedpe-like format will likeyly be the best way to do it
#   # how do I score a rearrangement?
#   # size scoring/filtering will likely be necessary
# }
# my $currClust; # MOVE ME WHEN DONE
# sub cdhitResults {
#   chomp($line);
#   if ($line =~ /^>(Cluster \d+)/) {
#     $currClust = $1;
#     $cluster{$1} = [];
#   } elsif ($line =~ /^\d+\t(.+)$/) {
#     my @f = split("\s", $1);
#     my $readName = substr($f[1], 1, length($f[1])-4)
#     if ($f[2] = "*") {
#
#     } elsif ($f[3] =~ ) {
#
#     }
#   }
# }
#
# cd-hit-est -A .95 -aS .95 -aL .95 -T 4 -c .99 -G 0 -mask N -i -o nr99
# cd-hit-est -A .95 -aS .95 -aL .95 -T 4 -c .95 -G 0 -mask N -i nr99 -o nr95
# cd-hit-est -A .95 -aS .95 -aL .95 -T 4 -c .85 -G 0 -mask N -i nr95 -o nr85 -n 6
# cd-hit-est -A .95 -aS .95 -aL .95 -T 4 -c .75 -G 0 -mask N -i nr85 -o nr75 -n 4

my %ambig = (
  "AG" => "R", "CT" => "Y", "GT" => "K", "AC" => "M", "CG" => "S", "AT" => "W", "CGT" => "B", "AGT" => "D", "ACT" => "H", "ACG" => "V"
);

# check that r1 r2 readnames can be matched
# if these reads overlap figure out how to select higher quality bases
# if they don't overlap use start/end positions to fill in Ns appropriately
