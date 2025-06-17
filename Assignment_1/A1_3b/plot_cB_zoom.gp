reset
set title 'Zeitverlauf von cB (bis 3s)'
set xlabel 't [s]'
set ylabel 'cB [mol/L]'
set xrange [0:3]
plot 'Ergebnisse.txt' using 1:3 with lines title 'cBlow', \
     '' using 1:6 with lines title 'cBhigh', \
     '' using 1:9 with lines title 'cBextr'
