=cut header

	running format:
		perl matching.pl prediction_complexes known_complexes
	function:
	    The function of this procedure is compute the the number of matching between prediction complexes and
	    known complexes.
	date:
	    07,21,2010
=cut

use warnings;
use strict;

my $firstPara = shift;			#prediction_complexes file
my $secondPara = shift;			#known complexes file
my $outputfile=shift;
open(my $prediction, $firstPara)||die("Could not open prediction_complexes file\n$!");
open(my $output, ">outputfile.txt")||die("Could not open output_matching.txt\n$!");
print $output "PComplex\tPSize\tKComplex\tKSize\tOverlap\t\tOS\n";
my $count = 0;	#the number of prediction complexes

while(<$prediction>){
	chomp $_;
	if($_ =~ m/Complex\s+(\d+)\s+(\d+)/){
		$count++;
		my $first = $1;
		my $num1 = $2;
		my @Pcomplex;
		for(my $i=1; $i<=$num1; $i++){
			my $line1 = <$prediction>;
			chomp $line1;
			push @Pcomplex, $line1;
		}
		my $Pcomplex = @Pcomplex;
		open(my $known, $secondPara)||die("Could not open known complexes file\n$!");
		while(<$known>){
			chomp $_;
			if($_ =~ m/Complex\s+(\d+)\s+(\d+)/){
				my $second = $1;
				my $num2 = $2;
				my @Kcomplex;
				for(my $j=1; $j<=$num2; $j++){
					my $line2 = <$known>;
					chomp $line2;
					push @Kcomplex, $line2;
				}
				my $Kcomplex = @Kcomplex;
				#compute the number of proteins between a prediction complex and a konwn complex which are same
				my $number = 0;
				for(my $m=0; $m<$num1; $m++){
					for(my $n=0; $n<$num2; $n++){
						if($Pcomplex[$m] eq $Kcomplex[$n]){
							$number++;
						}
					}
				}
				#compute the matching value, formular_1: M=(i**2)/(a*b)
				my $M1 = ($number**2)/($num1*$num2);
				if($M1 >= 0.1){
					print $output "Complex".$first."\t".$Pcomplex."\t"."Complex".$second."\t".$Kcomplex."\t".$number."\t".$M1."\n";
				}
			}
		}
		
		close($known);
	}
}
close($output);
close($prediction);

my @result;
push @result, "\n\nStatistical Information:\n";
push @result, "OS\t\tPc\t\tMPc\t\tMKc\t\t\tSn\t\t\t\t\t\t\tSp\t\t\t\t\t\t\tF\n";

for(my $threshold=0.1; $threshold<=1.0; $threshold+=0.1){
	my (@pc, @kc, @pcUniq, @kcUniq);
	open($output, "outputfile.txt")||die("Could not open output_matching.txt\n$!");
	while(<$output>){
		chomp $_;
		my @temp = split("\t", $_);
		if($temp[5] ge $threshold){
			push @pc, $temp[0];
			my %hashOne = ();
			@pcUniq = grep { ! $hashOne{$_} ++ } @pc;
			push @kc, $temp[2];
			my %hashTwo = ();
			@kcUniq = grep { ! $hashTwo{$_} ++ } @kc;
		}
	}
	my $numPc =  @pcUniq;
	my $numKc = @kcUniq;
	close($output);
	
	my $sn = $numPc/($numPc+(408-$numKc));
	my $sp = $numPc/($numPc+($count-$numPc));
	my $f;
	if($sn+$sp != 0){
		$f = 2*$sp*$sn/($sp+$sn); 
	}else{
		$f=0;
	}
	push @result, "$threshold\t\t$count\t\t$numPc\t\t$numKc\t\t$sn\t\t$sp\t\t$f\n";
}

open($output, ">>$outputfile")||die("Could not open output_matching.txt\n$!");
print $output $firstPara."---------------------------\n";
foreach(@result){
	print $output $_;
}
print $output "--------------------------------\n";
close($output);
