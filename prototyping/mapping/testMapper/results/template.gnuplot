
# sphere ----------------------------------------------------------------------------

set terminal png size 900, 600
set output 'sphere160x160_error.png'
splot '../meshes/sphere160x160_error.gnuplot' using 1:2:3:4 with points pt 10  palette

set terminal png size 900, 600
set output 'sphere80x80_error.png'
splot '../meshes/sphere80x80_error.gnuplot' using 1:2:3:4 with points pt 10  palette 


# diffBound ---------------------------------------------------------------------------

set terminal png size 900, 600
set output 'diffBoundMaster_error.png'
splot '../meshes/diffBoundMaster_error.gnuplot' using 1:2:3:4 with points pt 10  palette

set terminal png size 900, 600
set output 'diffBoundSlave_error.png'
splot '../meshes/diffBoundSlave_error.gnuplot' using 1:2:3:4 with points pt 10  palette 


# matchingDistortQuad ------------------------------------------------------------------

set terminal png size 900, 600
set output 'matchingDistortQuadMaster_error.png'
splot '../meshes/matchingDistortQuadMaster_error.gnuplot' using 1:2:3:4 with points pt 10  palette

set terminal png size 900, 600
set output 'matchingDistortQuadSlave_error.png'
splot '../meshes/matchingDistortQuadSlave_error.gnuplot' using 1:2:3:4 with points pt 10  palette 


# patchTest ------------------------------------------------------------------------------

set terminal png size 900, 600
set output 'patchTestMaster_error.png'
splot '../meshes/patchTestMaster_error.gnuplot' using 1:2:3:4 with points pt 10  palette

set terminal png size 900, 600
set output 'patchTestSlave_error.png'
splot '../meshes/patchTestSlave_error.gnuplot' using 1:2:3:4 with points pt 10  palette 


# sameBoundTri ----------------------------------------------------------------------------

set terminal png size 900, 600
set output 'sameBoundTriMaster_error.png'
splot '../meshes/sameBoundTriMaster_error.gnuplot' using 1:2:3:4 with points pt 10  palette

set terminal png size 900, 600
set output 'sameBoundTriSlave_error.png'
splot '../meshes/sameBoundTriSlave_error.gnuplot' using 1:2:3:4 with points pt 10  palette 





 