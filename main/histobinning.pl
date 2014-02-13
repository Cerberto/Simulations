#!/usr/bin/perl -w

###############################################################
#
# Script to bin data
#
# Run:
# ./binning.pl input.dat > output.dat
#
# input.dat has just 1 column, which are the value of the r.v..
#
# output.dat has 3 columns which are:
#	_ bin's lower value;
#	_ bin's upper value;
#	_ number of counts in the bin.
#
###############################################################


$step = 100.0;		# number of bins
$x_min = 0.5;		# lower value of the range
$x_max = 100.5;		# upper value of the range
$dx = ($x_max - $x_min)/$step;	# bin's width

@HISTO = (0) x $step;

$i = 0;
# count the data per bin and store in the HISTO array
while(<>){
	chomp;
	$i++;
	$j = int (($_ - $x_min)/$dx);
	$HISTO[$j]++;
}

# write data on output file (stdout)
for($i = 0; $i < $step; $i++){
	$x1 = $x_min + $i*$dx;
	$x2 = $x_min + ($i + 1)*$dx;
	printf "%e\t%e\t%d\n", $x1, $x2, $HISTO[$i];
}
