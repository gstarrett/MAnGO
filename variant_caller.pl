#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::DB::Sam;
use Getopt::Long;

my $infile = shift; # Input fastq files for alignment
my $ref = shift; # Reference fasta

my $aln;
# Path to bwa mem command
my $bwamem = "/Applications/bioinf_tools/bwa-0.7.5a/bwa mem";
my $samPath = "$infile.sam";
my $bamPath = "$infile.bam";
my $sortPath = "$infile.sort.bam";
my (@summary, @prefilter, @out, %readLenHash, %posHash, %readHash, %baseHash);

# system("$bwamem -B 2 $ref $infile > $samPath");
system("/Applications/bioinf_tools/samtools-0.1.19/samtools view -Sbt $ref.fai -o $bamPath $samPath");
system("/Applications/bioinf_tools/samtools-0.1.19/samtools sort -f $bamPath $sortPath");
system("/Applications/bioinf_tools/samtools-0.1.19/samtools index $sortPath");

#$sortPath = $infile;

### high level API
my $sam = Bio::DB::Sam->new(-bam => $sortPath, -fasta=> $ref);
my @targets = $sam->seq_ids;

my $snp_caller = sub {
			 my ($seqid,$pos,$p) = @_;
			 #print Dumper(@_);
			 my $refbase = uc($sam->segment($seqid,$pos,$pos)->dna);
			 my $trint = uc($sam->segment($seqid,$pos-1,$pos+1)->dna);
			 my $flank = uc($sam->segment($seqid,$pos-10,$pos+10)->dna);
			 next if $refbase =~ /[nN]/;
			 for my $pileup (@$p) {
					 my $b = $pileup->alignment;
					 my $inBool = $pileup->indel;
					 my $delBool = $pileup->is_refskip;
					 #print join("\t", $inBool, $delBool), "\n";
					 if ($inBool == 1 || $delBool == 1) {

					 } else {
						 my $qbase = uc(substr($b->qseq,$pileup->qpos,1));
						 next if $qbase =~ /[nN]/;
						 #print join("\t", $qbase, $refbase), "\n";
						 my $qscore = $b->qscore->[$pileup->qpos];
						 next unless $qscore >= 30;
						 my $qreadname = $b->query->name;
						 $baseHash{$qbase}++;
						 if($refbase ne $qbase) {
							 	$readHash{$qreadname}++;
								$posHash{$pos}++;
						 		push (@prefilter, join("\t",$qreadname,$seqid,$pos,$refbase,$qbase,$qscore,$trint,$flank));
					 	 }
				 	}
			 }
	 };

$sam->pileup($targets[0],$snp_caller);

open(SUMMARY, "> $infile.summary.txt");
# print SUMMARY "Number of reads used:\t$i\n";

my $posOutliers = calculate(\%posHash);
my $readOutliers = calculate(\%readHash);

my @posFlt;

print SUMMARY "===Unfiltered base counts===\n";
for my $baseKey (keys %baseHash) {
	print SUMMARY "$baseKey:\t$baseHash{$baseKey}\n";
}
print SUMMARY "===Filtered base counts===\n";
foreach (@prefilter) {
	if (scalar @$posOutliers == 0) {
		push(@posFlt, "$_\tPASS");
	} else {
		my $bool = 0;
		for my $posOutlier (@$posOutliers) {
			if ($bool == 0) {
				if ($_ =~ /$targets[0]\t$posOutlier\t/) {
					$bool = 1;
				} else {
					next;
				}
			} else {
				next;
			}
		}
		if ($bool == 0 ) {
			push(@posFlt, "$_\tPASS");
		}
		else {
			push(@posFlt, "$_\tOUTLIER");
		}
	}
}

foreach (@posFlt) {
	if (scalar @$readOutliers == 0) {
		push(@out, "$_\tPASS");
	} else {
		my $bool = 0;
		for my $readOutlier (@$readOutliers) {
			if ($bool == 0) {
				if ($_ =~ /^$readOutlier\t/) {
					$bool = 1;
				} else {
					next;
				}
			} else {
				next;
			}
		}
		if ($bool == 0) {
			push(@out, "$_\tPASS");
		} else {
			push(@out, "$_\tOUTLIER");
		}
	}
}

if (scalar @$readOutliers == 0) {
	for my $baseKey (keys %baseHash) {
		print SUMMARY "$baseKey:\t$baseHash{$baseKey}\n";
	}
} else {
	for my $readOutlier (@$readOutliers) {
		my $int = $readLenHash{$readOutlier};
		for my $baseKey (keys %{$int}) {
			$baseHash{$baseKey} = $baseHash{$baseKey} - ${$int}{$baseKey};
		}
	}
	for my $baseKey (keys %baseHash) {
		print SUMMARY "$baseKey:\t$baseHash{$baseKey}\n";
	}
}

open(OUT, "> $infile.counts.txt");

foreach (@out) {
	print OUT $_, "\n";
}


close(OUT);
close(SUMMARY);

#############
#SUBROUTINES#
#############

# Calculate if any read has unusually high (u+3sd) frequency of mutations and remove
# Calculate if any position has unusually high (u+3sd) frequency of mutations
sub calculate {
	my $hash = shift;
	#print Dumper(%{$hash});
	my $sum = 0;
	my $sumDev = 0;
	my @outliers;
	for my $key (keys %$hash) {
		$sum += $$hash{$key};
	}
	my $total = (scalar keys %$hash);
	my $mean = $sum/$total;
	for my $key (keys %$hash) {
		$sumDev += ($$hash{$key} - $mean) ** 2;
	}
	my $stdev = sqrt($sumDev/$total);
	for my $key (keys %$hash) {
		if ($$hash{$key} > ($mean + 3*$stdev)) {
			push(@outliers, $key);
		}
	}
	return(\@outliers);
}
