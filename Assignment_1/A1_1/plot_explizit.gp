reset
set title 'Explizites Euler-Verfahren'
set xlabel 't [s]'
set ylabel 'c [mol/L]'
plot 'Explizit_Ergebnisse.txt' using 1:2 with lines title 'cA', 'Explizit_Ergebnisse.txt' using 1:3 with lines title 'cB', 'Explizit_Ergebnisse.txt' using 1:4 with lines title 'cC'
