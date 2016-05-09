#!/usr/bin/gnuplot

radius = 6378100.0
sx(u,v) = radius*cos(u)*cos(v)
sy(u,v) = radius*cos(u)*sin(v)
sz(u,v) = radius*sin(u)

range = 10*radius

set parametric

set urange [-0.5*pi:0.5*pi]
set vrange [0:2*pi]

set xrange [-range:range]
set yrange [-range:range]
set zrange [-range:range]
set xlabel "x"
set ylabel "y"
set zlabel "z"

splot 'orbit.dat', sx(u,v),sy(u,v),sz(u,v)

pause -1

