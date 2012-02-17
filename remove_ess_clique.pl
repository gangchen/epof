use warnings;
use strict;

my $complexFilename = shift;
my %clusteredNodes;
open(my $complexFile, $complexFilename) or die $!;
while(my $line = <$complexFile>){
    next if($line =~ /Complex/);
    chomp($line);
    $clusteredNodes{$line} = 1;
}
close $complexFile;

open(my $netFile, "weight.txt") or die $!;
while(my $line = <$netFile>){
    chomp($line);
    my ($proteinA, $proteinB) = split "\t",$line;
    if(defined($clusteredNodes{$proteinA}) and 
       defined($clusteredNodes{$proteinB}))
    {
	print $line, "\n";
    }
}
close $complexFile;
