#!/usr/bin/perl -w

$i = 0;
$N = 12;

for ($i=0; $i < $N; $i+=1)
{
	system("./lennard-jones-P < input/LJ_P_input.in");
}

printf "\nEnd!\n\n";
