// Paramètres du maillage
meshEp = 0.5;
numLayers = 10;

// Points du carré
Point(1) = {-30.0, -10.0, 0.0, meshEp};
Point(2) = {10.0, -10.0, 0.0, meshEp};
Point(3) = {10.0, 10.0, 0.0, meshEp};
Point(4) = {-30.0, 10.0, 0.0, meshEp};

// Lignes qui relient les points
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Carré interne
Point(5) = {-21.5, -1.5, 0.0, meshEp};
Point(6) = {-18.5, -1.5, 0.0, meshEp};
Point(7) = {-18.5, 1.5, 0.0, meshEp};
Point(8) = {-21.5, 1.5, 0.0, meshEp};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// Ligne pour déterminer la surface du carré
Line Loop(9) = {1,2,3,4};
Line Loop(10) = {5,6,7,8};

// Surface du carre
Plane Surface(1) = {9,10};

// Géométrie 3D
// --> Pas nécéssaire: Extrude {0,0,0.1}{Surface{1};Layers{numLayers};Recombine;}

