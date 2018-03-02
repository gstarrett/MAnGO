#!/usr/bin/perl -w
use strict;

my (%clstrHash1, %clstrHash2, %clstrHash3, $clstrName);
my $numLevels = 3; # 999, 99, 95
my $infile1 = shift;
my $infile2 = shift;
my $infile3 = shift;
open (CLSTR, "< $infile1");
while (my $line = <CLSTR>) {
  chomp($line);
  my $readName;
  if ($line =~ />([a-z0-9]+)\.\.\./) {
    $readName = $1;
    # print "Found $1!\n";
    if ($line =~ /\*/) {
      $clstrName = $readName;
    }
    push(@{$clstrHash1{$clstrName}}, $readName);
  }
}
close(CLSTR);

open (CLSTR, "< $infile2");
while (my $line = <CLSTR>) {
  chomp($line);
  my $readName;
  if ($line =~ />([a-z0-9]+)\.\.\./) {
    $readName = $1;
    # print "Found $1!\n";
    if ($line =~ /\*/) {
      $clstrName = $readName;
    }
    push(@{$clstrHash2{$clstrName}}, $readName);
  }
}
close(CLSTR);

open (CLSTR, "< $infile3");
while (my $line = <CLSTR>) {
  chomp($line);
  my $readName;
  if ($line =~ />([a-z0-9]+)\.\.\./) {
    $readName = $1;
    # print "Found $1!\n";
    if ($line =~ /\*/) {
      $clstrName = $readName;
    }
    push(@{$clstrHash3{$clstrName}}, $readName);
  }
}
close(CLSTR);

# collect seqnames under clusters, provide representative seqname to clusters
#repeat for next level

for my $key3 (keys %clstrHash3) {
  for my $key2 (@{$clstrHash3{$key3}}) {
    for my $key1 (@{$clstrHash2{$key2}}) {
      for my $name (@{$clstrHash1{$key1}}) {
        print join("\t", $key3, $key2, $key1, $name), "\n";
      }
    }
  }
}

# data.wide.0.1.alpha <- fisher.alpha(round(data.wide.0.1[,5:ncol(data.wide.0.1)]*12000), MARGIN = 2)
# data.wide.0.2.alpha <- fisher.alpha(round(data.wide.0.2[,5:ncol(data.wide.0.2)]*12000), MARGIN = 2)
# data.wide.0.3.alpha <- fisher.alpha(round(data.wide.0.3[,5:ncol(data.wide.0.3)]*12000), MARGIN = 2)
# data.wide.clstr1.0.1.alpha <- fisher.alpha(round(data.wide.clstr1.0.1[,4:ncol(data.wide.clstr1.0.1)]*12000), MARGIN = 2)
# data.wide.clstr1.0.2.alpha <- fisher.alpha(round(data.wide.clstr1.0.2[,4:ncol(data.wide.clstr1.0.2)]*12000), MARGIN = 2)
# data.wide.clstr1.0.3.alpha <- fisher.alpha(round(data.wide.clstr1.0.3[,4:ncol(data.wide.clstr1.0.3)]*12000), MARGIN = 2)
# data.wide.clstr2.0.1.alpha <- fisher.alpha(round(data.wide.clstr2.0.1[,3:ncol(data.wide.clstr2.0.1)]*12000), MARGIN = 2)
# data.wide.clstr2.0.2.alpha <- fisher.alpha(round(data.wide.clstr2.0.2[,3:ncol(data.wide.clstr2.0.2)]*12000), MARGIN = 2)
# data.wide.clstr2.0.3.alpha <- fisher.alpha(round(data.wide.clstr2.0.3[,3:ncol(data.wide.clstr2.0.3)]*12000), MARGIN = 2)
# data.wide.clstr2.all.alpha <- rbind(data.wide.clstr2.0.1.alpha, data.wide.clstr2.0.2.alpha, data.wide.clstr2.0.3.alpha)
# data.wide.clstr1.all.alpha <- rbind(data.wide.clstr1.0.1.alpha, data.wide.clstr1.0.2.alpha, data.wide.clstr1.0.3.alpha)
# data.wide.all.alpha <- rbind(data.wide.0.1.alpha, data.wide.0.2.alpha, data.wide.0.3.alpha)

# ggplot(melt(data.wide.clstr2.all.alpha), aes(x= Var2, y=value, shape=c("A","B","C","A","B","C","A","B","C","A","B","C","A","B","C","A","B","C","A","B","C","A","B","C","A","B","C","A","B","C"), color=c("A3B","A3B","A3B","Control","Control","Control","A3B","A3B","A3B","Control","Control","Control","A3B","A3B","A3B","Control","Control","Control","A3B","A3B","A3B","Control","Control","Control","A3B","A3B","A3B","Control","Control","Control"))) + geom_line() + geom_point() + theme(legend.position = "none")
# plot.con.3 <- ggplot(data.long[which(data.long$exp == "control" & data.long$rep == 3),], aes(week, ID, fill = freq)) + geom_raster(interpolate=TRUE) + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),legend.position="none") + scale_fill_gradientn(colours=c("white","goldenrod1","orangered","darkred","black"), limits=c(0,0.1), values = c(0,0.002,0.012,0.06,1))
