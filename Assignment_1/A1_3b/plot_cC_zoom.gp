reset
set title 'Zeitverlauf von cC (bis 3s)'
set xlabel 't [s]'
set ylabel 'cC [mol/L]'
set xrange [0:3]
plot 'Ergebnisse.txt' using 1:4 with lines title 'cClow', \
     '' using 1:7 with lines title 'cChigh', \
     '' using 1:10 with lines title 'cCextr'
