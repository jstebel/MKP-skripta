#!/usr/bin/gnuplot

set term pdf enhanced
set output 'discontinuous_h1.pdf'
set pm3d
set isosample 500
set view 60,300
set ztics 0.5
splot [-1:1][-1:1] (x*x+y*y>1)?1/0:(x*x+y*y)**(0.25)*sin(atan2(x,y)/2) with pm3d not