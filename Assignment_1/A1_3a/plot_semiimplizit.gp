reset
set title 'Semi-Implizites Euler + Richardson'
set xlabel 't [s]'
set ylabel 'c [mol/L]'
plot 'SemiImplizit_Richardson.txt' using 1:2 with lines title 'cA', 'SemiImplizit_Richardson.txt' using 1:3 with lines title 'cB', 'SemiImplizit_Richardson.txt' using 1:4 with lines title 'cC'
