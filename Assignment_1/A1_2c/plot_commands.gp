reset
set terminal qt
set xlabel 't [s]'
set ylabel 'Konzentration [mol/L]'
set logscale y
set yrange [1e-8:1.2]
set title 'Explizit Euler, h=2e-6'
plot 'Ergebnisse_Explizit_h=2e-6.txt' using 1:2 with lines title 'cA', 'Ergebnisse_Explizit_h=2e-6.txt' using 1:3 with lines title 'cB', 'Ergebnisse_Explizit_h=2e-6.txt' using 1:4 with lines title 'cC'
pause mouse close
set title 'Explizit Euler, h=1e-6'
plot 'Ergebnisse_Explizit_h=1e-6.txt' using 1:2 with lines title 'cA', 'Ergebnisse_Explizit_h=1e-6.txt' using 1:3 with lines title 'cB', 'Ergebnisse_Explizit_h=1e-6.txt' using 1:4 with lines title 'cC'
pause mouse close
set title 'Explizit Euler, h=5e-7'
plot 'Ergebnisse_Explizit_h=5e-7.txt' using 1:2 with lines title 'cA', 'Ergebnisse_Explizit_h=5e-7.txt' using 1:3 with lines title 'cB', 'Ergebnisse_Explizit_h=5e-7.txt' using 1:4 with lines title 'cC'
pause mouse close
set title 'Implizit Euler, h=2e-6'
plot 'Ergebnisse_Implizit_h=2e-6.txt' using 1:2 with lines title 'cA', 'Ergebnisse_Implizit_h=2e-6.txt' using 1:3 with lines title 'cB', 'Ergebnisse_Implizit_h=2e-6.txt' using 1:4 with lines title 'cC'
pause mouse close
set title 'Implizit Euler, h=1e-6'
plot 'Ergebnisse_Implizit_h=1e-6.txt' using 1:2 with lines title 'cA', 'Ergebnisse_Implizit_h=1e-6.txt' using 1:3 with lines title 'cB', 'Ergebnisse_Implizit_h=1e-6.txt' using 1:4 with lines title 'cC'
pause mouse close
set title 'Implizit Euler, h=5e-7'
plot 'Ergebnisse_Implizit_h=5e-7.txt' using 1:2 with lines title 'cA', 'Ergebnisse_Implizit_h=5e-7.txt' using 1:3 with lines title 'cB', 'Ergebnisse_Implizit_h=5e-7.txt' using 1:4 with lines title 'cC'
pause mouse close
