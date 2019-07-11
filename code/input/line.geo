lc = 0.004;
//+
Point(1) = {0.0, 0.0, 0, lc};
//+
Point(2) = {1, 0, 0, lc};
//+
Point(3) = {1, 0.005, 0, lc};
Point(4) = {0, 0.005, 0, lc};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Physical Curve("walls") = {1, 3};
//+
Physical Curve("inout") = {4, 2};
