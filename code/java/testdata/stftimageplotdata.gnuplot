set title "STFT and maximiser"
unset key
set tic scale 0

# Color runs from white to green
set palette rgbformula -7,2,-7
set cbrange [0:5]
set cblabel "Magnitude"
unset cbtics

set xrange [0:255]
set yrange [-0.5:0.5]

set view map
splot 'stftimageplotdata' matrix with image