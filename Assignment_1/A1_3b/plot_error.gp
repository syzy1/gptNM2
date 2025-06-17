reset
set title 'Fehler gegen low/extr und high/extr'
set xlabel 't [s]'
set ylabel 'Error'
plot 'Ergebnisse.txt' using 1:11 with lines title 'Errlow', \
     '' using 1:12 with lines title 'Errhigh'
