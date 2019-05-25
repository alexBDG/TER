lc =0.3;

// generation du rectangle
Point(1) = {-30,-10,0,lc};
Point(2) = {-30,10,0,lc};
Point(3) = {10,10,0,lc};
Point(4) = {10,-10,0,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// generation de l’ellipse
Point(5) = {-20,0,0,lc};
Point(6) = {-17,0,0,lc};
Point(7) = {-20,1.5,0,lc};
Point(8) = {-23,0,0,lc};
Point(9) = {-20,-1.5,0,lc};

Ellipse(5) = {6,5,6,7};
Ellipse(6) = {7,5,8,8};
Ellipse(7) = {8,5,8,9};
Ellipse(8) = {9,5,6,6};

// generation des contours
Line Loop(9) = {1,2,3,4};
Line Loop(10) = {5,6,7,8};

// generation de la surface
Plane Surface(1) = {1,2};
