Here is some explanation of the mortarmapper related work.

1. Mortar mapper 
	* Properties: 
		- consistent and conservative at the same time (only when with the same integration domain, e.g. a plane)
		- mapping on surfaces with curvature, there is problem both on the boundary and the inner domain
		- if a quadrilateral is distorted, than the Gauss quadrature may be not exact, Quadmortarmapper shows error
	* Test cases are designed for showing the above properties, they can be run by scripts inside folder "startScripts"
		- plane.sh                does mapping on a plane, where both consistent and conservative mapping are obtained
		- sphere.sh               does mapping on a sphere, the error can be observed on the boundary
		- diffBound.sh            shows the error due to nonmatching boundary
		- sameBoundTri.sh         shows the error due to different resolution of the same curvature
		- matchingDistortQuad.sh  shows the error due to error of Gauss quadrature when the quadrilateral is distorted
		- patchTest.sh            shows does mapping on a simple cylinder which is a classical patch test
		
2. Futrual Work:
	* Simplify code in GiDParser, redesign the API, add comments
	* Output pressure instead of forces (or both can be output)
