#!/usr/bin/perl -w

# Script for repeating several times the computation of linpot_lattice
# for several values of the variational parameters

$g = 0;
$gmin	= 0.05;
$gmax	= 5.8;
$N 		= 50;
$dg		= ($gmax - $gmin)/$N;

printf "Minimum -- maximum variational parameter : $gmin -- $gmax\n";

system ("rm", "-f", "lp_output/expectationvalues_bf.dat");
system ("touch", "lp_output/expectationvalues_bf.dat");

for ($g=$gmin; $g < $gmax; $g+=$dg)
{
	# printf "Variational parameter = $g\n";
	system("./linpot_lattice $g < input/LP_input_bf.in");
}

printf "\nEnd!\n\n";