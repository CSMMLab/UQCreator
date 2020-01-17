// scaling factors. MODIFY LAMBDA for individual problem
lambda = 92.6*1e-6*100;
//lambda = 9.26*1e-6*100;
//lambda = 15.53*1e-6*100;
sigma = 1.0/lambda;
scale = 10.0;

// resolution
//lc = 0.001*sigma;
//lcInner = 0.0001*sigma;
lc = 0.001*sigma*scale;
lcInner = 0.0001*sigma*scale;
//lc = 0.01*sigma;
//lcInner = 0.001*sigma;

// domain bounds
a = -0.02*sigma*scale;
b = 0.02*sigma*scale;
//a = -0.2*sigma;
//b = 0.2*sigma;
aInner = (-2.5e-4 * 0.5*sigma)*scale;
bInner = (2.5e-4 * 0.5*sigma)*scale;

//+
Point(1) = {b, a, 0, lc};
//+
Point(2) = {b, b, 0, lc};
//+
Point(3) = {a, b, 0, lc};
//+
Point(4) = {a, a, 0, lc};
//+
Point(5) = {bInner, aInner, 0, lcInner};
//+
Point(6) = {bInner, bInner, 0, lcInner};
//+
Point(7) = {aInner, bInner, 0, lcInner};
//+
Point(8) = {aInner, aInner, 0, lcInner};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {7, 6};
//+
Line(6) = {6, 5};
//+
Line(7) = {5, 8};
//+
Line(8) = {8, 7};

//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Curve Loop(2) = {5, 6, 7, 8};
//+
Plane Surface(1) = {1, 2};
//+
Plane Surface(2) = {2};

//+
Physical Curve("OuterBoundary") = {1, 4, 3, 2};
