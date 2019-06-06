
lc = 0.025;
a = -1.5;
b = 1.5;
 //+
Point(1) = {b, a, 0, lc};
//+
Point(2) = {b, b, 0, lc};
//+
Point(3) = {a, b, 0, lc};
//+
Point(4) = {a, a, 0, lc};
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
Physical Curve("OuterBoundary") = {1, 4, 3, 2};
