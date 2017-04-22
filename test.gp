set encoding utf8
set key outside Left
set bmargin 5
set tmargin 6
set style data lines
set tics in
set ticslevel 0.5
set xlabel  "X-ray energy in eV"
set format y  '%5.1fe'
set title " Anomalous scattering factors "
set xrange  [9000:14400]
set offset 0,0,1.0,0
set xtics nomirror
set link x via 12398./x inverse 12398./x

set x2label  "X-ray wavelength in Å"
set x2tics 0.1  format "%.1f Å" nomirror

Brdata = "< GET http://www.bmsc.washington.edu/scatter/data/Br.dat"
Tadata = "< GET http://www.bmsc.washington.edu/scatter/data/Ta.dat"
pause -1 "Hit any key to continue"
