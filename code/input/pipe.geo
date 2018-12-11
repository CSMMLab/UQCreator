lcC = 2.0; // resolution farfield
lcF = 0.015;

// inner rectangle points
Point(1)={0, 0, 0, lcF};
Point(2)={1, 0, 0, lcF};
Point(3)={2, 1, 0, lcF};
Point(4)={2, 2, 0, lcF};
Point(5)={2, 3, 0, lcF};

// inner rectangle
//Line(1)={1,2};
//Line(2)={2,3};
//Line(3)={3,4};
//Line(4)={4,5};
//+
Point(6) = {1.8, 0.2, 0, lcF};
//+
Point(7) = {0, 0.5, -0, lcF};
//+
Point(8) = {1, 0.5, 0, lcF};
//+
Point(9) = {1.4, 0.6, 0, lcF};
//+
Point(10) = {1.5, 1, 0, lcF};
//+
Point(11) = {1.5, 2, 0, lcF};
//+
Point(12) = {1.5, 3, 0, lcF};
//+
Spline(1) = {1, 2, 6, 3, 4, 5};
//+
Spline(2) = {7, 8, 9, 10, 11, 12};
//+
Line(3) = {7, 1};
//+
Line(4) = {12, 5};
//+
Curve Loop(1) = {3, 1, -4, -2};
//+
Plane Surface(1) = {1};
//+
Physical Curve("inlet") = {4};
//+
Physical Curve("outlet") = {3};
//+
Physical Curve("walls") = {2, 1};
