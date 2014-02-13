
reset
set terminal epslatex size 12.5cm,7.5cm color colortext
set output 'histo_DeltaE.tex'

set format '$%g$'
unset key
set ylabel 'Frequency'
set xlabel 'Site'
set xrange [0:15]
set yrange [0:50000]
#set y2label 'Cumulative frequency'
#set y2tics 0,2e3,1e4	# parte da 0, salta di 2000 e arriva a 10^4
set ytics nomirror
set boxwidth 0.8 relative	# larghezza box relativa alla larghezza del bin
set style fill solid 0.4	# riempimento con trasparenza 60%


a = 45000.0
g = 0.4
f(x) = a*exp(-g*x)
set label 1 "$a = %.4g$", a at (2/g),(a-a/3)
set label 2 "$g = %.3g$", g at (2/g),(a-a/2)

set style line 1 lt 1 lc rgb '#0060ad' lw 3 pt 1

fit f(x) 'distribution_binned.dat' using (($1+$2)/2):3 via a,g
plot 'distribution_binned.dat' using (($1+$2)/2):3 with impulses, \
	f(x) with lines ls 1

#
#	Npoints=3000.0
#
#	bw = 0.001
#	bin(x,width)=width*floor(x/width)
#
#	plot 'deltaE_results.dat' using (bin($1,bw)):(1/Npoints) smooth frequency with impulses lc rgb '#dd181f'

