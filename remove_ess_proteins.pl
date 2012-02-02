# remove essential proteins from the network

my %essProteins;
open(my $essProteinFile, "key.txt") or die $!;
while(my $line = <$essProteinFile>){
    chomp($line);
    $essProteins{$line} = 1;
}
close($essProteinFile);

open(my $networkFile, "SC_net_unweight.txt") or die $!;
while(my $line = <$networkFile>){
    chomp($line);
    my ($proteinA, $proteinB, $weight) = split "\t", $line;
    if(defined($essProteins{$proteinA}) and 
       defined($essProteins{$proteinB})){
	print $line, "\n";
    }
}
close $networkFile;
