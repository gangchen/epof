

my %ess;
open(my $essFile, "key.txt") or die $!;
while(my $line = <$essFile>){
    chomp($line);
    $ess{$line} = 1;
}
close $essFile;

open(my $complexFile, "complex_408.txt") or die $!;
my $complexId;
while(my $line = <$complexFile>){
    chomp($line);
    my $proteinNum;
    if($line =~ /Complex\t(.+)\t(.+)/){
	$complexId = $1;
	$proteinNum = $2;
	$flag = 0;
	for(1 .. $proteinNum){
	    $line = <$complexFile>;
	    chomp($line);
	    if(defined($ess{$line})){
		$flag = 1;
	    }
	}
	if($flag == 1){
	    print $complexId, "\t", $proteinNum,
"\n";
	}

    }else{
	print "ERROR\n";
    }
}
close $complexFile;
