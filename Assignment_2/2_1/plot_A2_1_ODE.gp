reset
set title 'Loesung ODE: Semi-implizites Euler + Extrapolation'
set xlabel 't [min]'
set ylabel 'Nm, Np, Vr'
plot 'Loesung_A2_1_ODE.txt' using 1:2 with lines title 'Nm', \
     'Loesung_A2_1_ODE.txt' using 1:3 with lines title 'Np', \
     'Loesung_A2_1_ODE.txt' using 1:4 with lines title 'Vr'
