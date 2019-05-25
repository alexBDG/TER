// Paramètres du maillage
meshEp = 0.18;//0.314;
numLayers = 10;

// Points du carré
Point(1) = {-10, -10, 0.0, meshEp};
Point(2) = {-10, 10, 0.0, meshEp};
Point(3) = {10, 10, 0.0, meshEp};
Point(4) = {10, -10, 0.0, meshEp};

// Lignes qui relient les points
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Ligne pour déterminer la surface du carré
Line Loop(5) = {1,2,3,4};

// Surface du carre
Plane Surface(1) = {5};

// Géométrie 3D
// --> Pas nécéssaire: Extrude {0,0,0.1}{Surface{1};Layers{numLayers};Recombine;}
