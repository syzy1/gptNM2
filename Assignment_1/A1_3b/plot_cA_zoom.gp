reset
set title 'Zeitverlauf von cA (bis 3s)'
set xlabel 't [s]'
set ylabel 'cA [mol/L]'
set xrange [0:3]
plot 'Ergebnisse.txt' using 1:2 with lines title 'cAlow', \
     '' using 1:5 with lines title 'cAhigh', \
     '' using 1:8 with lines title 'cAextr'
