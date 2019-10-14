#
# Various ways to create a 2D heat map from ascii data
#

set terminal pngcairo size 1920,1080
set output "plot.png"

set title "Heat Map for the Density of States"
unset key
set tic scale 0

# Color runs from white to green
#set palette rgbformula -7,2,-7
#set palette rgbformula -7,2,-7
set palette defined ( 0 '#FFFFFF',\
    	    	      1 '#FFFFFF',\
		      2 '#FDD0A2',\
		      3 '#FDAE6B',\
		      4 '#FD8D3C',\
		      5 '#F16913',\
		      6 '#D94801',\
		      7 '#8C2D04' )

# set cbrange [0:5]
set cblabel "DoS"
unset cbtics

set xrange [-0.025:9.3]
set yrange [-1.0:0.5]

set view map
set pm3d interpolate 1,1
splot "logder_scan_s" using 1:2:(($3)**(0.125)) with image

