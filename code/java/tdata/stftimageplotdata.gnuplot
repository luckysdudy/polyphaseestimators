set title "STFT and maximiser"
set xrange [0:255]
set yrange [-0.5:0.5]
set terminal pdf
set output "stftimage.pdf"

plot 'stftimageplotdata' binary matrix with image, \
"stftmax" with lines