# plot_A2_3.gnuplot
# Grafische Auswertung für Aufgabe 2.3
set terminal pngcairo size 1200,900 enhanced font 'Arial,12'
set output 'A2_3_multiplot.png'
set xlabel 'Zeit t [min]'
set ylabel 'Werte'
set grid
set multiplot layout 2,2 title 'Aufgabe 2.3 – Modell 1 vs. Modell 2'

set title 'Monomermenge N_m'
plot \
    'Model1_A2_3.txt' using 1:2 with lines lw 2 title 'Modell 1', \
    'Model2_A2_3.txt' using 1:2 with lines lw 2 title 'Modell 2'

set title 'Polymermenge N_p'
plot \
    'Model1_A2_3.txt' using 1:3 with lines lw 2 title 'Modell 1', \
    'Model2_A2_3.txt' using 1:3 with lines lw 2 title 'Modell 2'

set title 'Reaktorvolumen V_r [m^3]'
plot \
    'Model1_A2_3.txt' using 1:4 with lines lw 2 title 'Modell 1', \
    'Model2_A2_3.txt' using 1:4 with lines lw 2 title 'Modell 2'

set title 'Zulaufstrom \\dot{V}'
plot \
    'Model1_A2_3.txt' using 1:5 with lines lw 2 title 'Modell 1', \
    'Model2_A2_3.txt' using 1:5 with lines lw 2 title 'Modell 2'

unset multiplot
unset output
