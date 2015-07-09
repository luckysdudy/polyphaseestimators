set title "STFT and maximiser"

set terminal pdf
set output "stftimage.pdf"

set xrange [0:63]
set yrange [-0.5:0.5]
plot 'm2nowrappingdata' binary matrix with image, \
"m2nowrappingmax" with lines

set xrange [0:255]
plot 'm3case1data' binary matrix with image, \
"m3case1max" with lines

set xrange [0:198]
plot 'm3case2data' binary matrix with image, \
"m3case2max" with lines

