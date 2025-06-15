set term qt 0 title 'Konzentrationen (semi-implizit)'
plot 'Fehler_Ergebnis.txt' using 1:2 with lines title 'cA', \
     '' using 1:3 with lines title 'cB', \
     '' using 1:4 with lines title 'cC'

set term qt 1 title 'Fehler E_i = ||Y_high−Y_low||'
plot 'Fehler_Ergebnis.txt' using 1:5 with lines title 'Fehler'
