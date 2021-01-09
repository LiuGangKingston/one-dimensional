# one-dimensional


Here is about the data files regarding the code one_dim.f90 and drawing graphs with the 
Gnuplot for the one-dimensional simplified model system in the paper:
Liu, G. A new equation for period vectors of crystals under external stress and temperature
in statistical physics: mechanical equilibrium condition and equation of state.
Eur. Phys. J. Plus 136, 48 (2021). https://doi.org/10.1140/epjp/s13360-020-01010-6 
January 08th, 2021.


The file in.dat is the only one input data file. It contains five lines of input data only, 
as follows. 

The first line is an integer, the total number of temperatures for the Constant external 
temperature cases, although the variable to get it is called total_external_force_lines in 
the code.

The second line is real numbers, all the temperature values for the Constant external 
temperature cases. 

The third line is one real number, the starting/minimum value of the external force, also 
for the Constant external temperature cases. 

The fourth line is an integer, the total number of the external forces for cases of Constant 
external force.

The last line is real numbers, all the external force values for cases of Constant external 
force.


Supposing all generated data files are in the directory "c:\one_dim_data\", curves in Figures 2 
through 7 in the paper can be drawn with the following commands inside Gnuplot. In the commands, 
following the "plot" is the data file to be read. Following "using", the two integers in the form 
of I:J are the two column numbers of the data to be used as x-component and y-component of the 
points in the graph respectively. All the following draws should be exactly the same as the 
corresponding Figures in the paper.


(Start Gnuplot for Figure 2)

plot 'c:\one_dim_data\single.forces.dat' using 1:2 title ' ' with lines linecolor "#FF0000", \
     'c:\one_dim_data\left.on.right.forces.dat' using 1:2 title ' ' with lines linecolor "#0000FF" , \
     'c:\one_dim_data\single.forces.dat' using 1:3 title ' ' with lines linecolor "#000000"


(Re-start Gnuplot for Figure 2_in)

set xrange [2.7:3.5]     

set tics font ", 20"

unset ytics

plot 'c:\one_dim_data\single.forces.dat' using 1:2 notitle with lines linecolor "#FF0000", \
     'c:\one_dim_data\left.on.right.forces.dat' using 1:2 notitle with lines linecolor "#0000FF" , \
     'c:\one_dim_data\single.forces.dat' using 1:3 notitle with lines linecolor "#000000"


(Re-start Gnuplot for Figure 3)

set yrange [2.6:3.2]  

plot 'c:\one_dim_data\period.vs.force.A.dat' using 2:3 title ' ' with lines linecolor "#FF0000", \
     'c:\one_dim_data\period.vs.force.B.dat' using 2:3 title ' ' with lines linecolor "#0000FF", \
     'c:\one_dim_data\period.vs.force.C.dat' using 2:3 title ' ' with lines linecolor "#00FF00"



(Re-start Gnuplot for Figure 4)

plot 'c:\one_dim_data\work.and.heat.A.dat' using 3:4 notitle  with lines linecolor "#FF0000", \
     'c:\one_dim_data\work.and.heat.B.dat' using 3:4 notitle  with lines linecolor "#0000FF", \
     'c:\one_dim_data\work.and.heat.C.dat' using 3:4 notitle  with lines linecolor "#00FF00"


(Re-start Gnuplot for Figure 5)
plot 'c:\one_dim_data\work.and.heat.A.dat' using 3:8 notitle  with lines linecolor "#FF0000", \
     'c:\one_dim_data\work.and.heat.B.dat' using 3:8 notitle  with lines linecolor "#0000FF", \
     'c:\one_dim_data\work.and.heat.C.dat' using 3:8 notitle  with lines linecolor "#00FF00"


(Re-start Gnuplot for Figure 6)

plot 'c:\one_dim_data\period.vs.temperature.A.dat' using 1:3 notitle  with lines linecolor "#FF0000", \
     'c:\one_dim_data\period.vs.temperature.B.dat' using 1:3 notitle  with lines linecolor "#0000FF", \
     'c:\one_dim_data\period.vs.temperature.C.dat' using 1:3 notitle  with lines linecolor "#00FF00"


(Re-start Gnuplot for Figure 7)

plot 'c:\one_dim_data\period.vs.temperature.A.dat' using 1:5 notitle  with lines linecolor "#FF0000", \
     'c:\one_dim_data\period.vs.temperature.B.dat' using 1:5 notitle  with lines linecolor "#0000FF", \
     'c:\one_dim_data\period.vs.temperature.C.dat' using 1:5 notitle  with lines linecolor "#00FF00"

