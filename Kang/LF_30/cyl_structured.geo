//Inputs by Felipe Huerta (2018)
// Imperial College London
radius = 0.0648;
height = 0.7581;
LF = 0.3;

// Grid size parameters

gridsize = radius*0.01;
dr = 1e-3;
dz = 1e-3;

Point(1) = {0,0, height*LF/2};
Point(2) = {0,0,-height*LF/2};
Point(3) = {radius,0, -height*LF/2};
Point(4) = {radius,0, height*LF/2};

Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,4};
Line(8) = {4,1};

// Calculate number of nodes in radial and axial directions
n_r = Round((radius)/dr);
n_z = Round((height*LF)/dz);

// Make non-uniform mesh in R

ref = 0.75;

Transfinite Line{6} = n_r Using Bump ref;
Transfinite Line{8} = n_r Using Bump ref;

Transfinite Line{7} = n_z Using Bump ref;
Transfinite Line{5} = n_z Using Bump ref;

Line Loop(9) = {5:8};
Plane Surface(10) = {9};
Transfinite Surface{10} = {1,2,3,4};

Rotate {{0,0,1},{0,0,0},2.5*Pi/180.0}
{
	Surface{10};
}

Recombine Surface{10};

new_entities[] = Extrude {{0,0,1},{0,0,0},-5*Pi/180.0}
{
	Surface{10};
	Layers{1};
	Recombine;
};

Physical Surface("wedge0") = {10};
Physical Surface("wedge1") = {new_entities[0]};
Physical Surface("bottom") = {new_entities[2]};
Physical Surface("tank_wall") = {new_entities[3]};
Physical Surface("interphase") = {new_entities[4]};

Physical Volume(1000) = {new_entities[1]};
