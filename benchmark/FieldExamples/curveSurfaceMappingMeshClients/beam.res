GiD Post Results File 1.0

GaussPoints "GP_triangles" ElemType Triangle "mesh_triangles"
  Number Of Gauss Points: 1
  Natural Coordinates: internal
end gausspoints

GaussPoints "GP_quads" ElemType Quadrilateral "mesh_quads"
  Number Of Gauss Points: 1
  Natural Coordinates: internal
end gausspoints

Result	"dispRot_disp"	"EMPIRE_CoSimulation"	1	vector	OnNodes
Values
	1	0	0	0
	2	1.03553	1.03553	-1.03553
	3	0	3.53553	-3.53553
	4	-3.96447	6.03553	-6.03553
	5	-10	7.07107	-7.07107
End Values

Result	"dispRot_rot"	"EMPIRE_CoSimulation"	1	vector	OnNodes
Values
	1	0	0	0
	2	0	0.55536	0.55536
	3	0	1.11072	1.11072
	4	0	1.66608	1.66608
	5	0	2.22144	2.22144
End Values

