=head1 
Script for sorting the nodes which is not included in any complex identified by the original EPOF algorithm
=cut


#Two parameters are needed: the file of the network, and the result file from EPOF
my $networkFilename, $complexFilename;

$complexFilename = shift;
$networkFilename = shift;

my %clusteredNodes;
open(my $complexFile, $complexFilename) or die $!;
while(my $line = <$complexFile>){
    next if($line =~ m/Complex/);
    chomp($line);
    $clusteredNodes{$line} = 1;
}
close $complexFile;

#print scalar keys %clusteredNodes;

my %degree;
open(my $networkFile, $networkFilename) or die $!;
while(my $line = <$networkFile>){
    chomp($line);
    my ($proteinA, $proteinB, $weight) = split "\t",$line;
    $degree{$proteinA}++;
    $degree{$proteinB}++;
}
close $networkFile;

my @sortedKey = sort {$degree{$b} <=> $degree{$a}} keys %degree;

for(@sortedKey){
    print $_, "\n" if(!defined($clusteredNodes{$_}));
}


