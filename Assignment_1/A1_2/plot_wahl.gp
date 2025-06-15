reset
set title 'Euler-Verfahren (Implizit)'
set xlabel 't [s]'
set ylabel 'c [mol/L]'
plot 'Wahl_Implizit.txt' using 1:2 with lines title 'cA', 'Wahl_Implizit.txt' using 1:3 with lines title 'cB', 'Wahl_Implizit.txt' using 1:4 with lines title 'cC'
