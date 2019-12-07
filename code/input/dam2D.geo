// Gmsh project created on Fri Dec  7 09:45:02 2018

// coarseness definition
lcCoarse = 0.03; // 0.05
lcFine = 0.005;

counter = 1;

// Geometry settings
posDamX = 0.5;
valueIntermediate = 0.2;
xIntermediate = 0.3;
posBreakYIntermLow = 0+valueIntermediate;
posBreakYLow = 0.35;
posBreakYUp = posBreakYLow+0.3;
posBreakYIntermUp = 1.0-valueIntermediate;
thicknessWall = 0.03;

// setup outer frame
Point(counter)={0.0,0.0,0, lcCoarse };
counter = counter+1;

Point(counter)={0.0,1.0,0, lcCoarse };
counter = counter+1;

Point(counter)={1.0,1.0,0, lcCoarse };
counter = counter+1;

Point(counter)={1.0,0.0,0, lcCoarse };
counter = counter+1;

// setup dam wall

// below break

Point(counter)={posDamX,0.0,0, lcCoarse };
counter = counter+1;

Point(counter)={posDamX,posBreakYLow,0, lcFine };
counter = counter+1;

Point(counter)={posDamX+thicknessWall,0.0,0, lcCoarse };
counter = counter+1;

Point(counter)={posDamX+thicknessWall,posBreakYIntermLow,0, lcFine };
counter = counter+1;

Point(counter)={posDamX+thicknessWall,posBreakYLow,0, lcFine };
counter = counter+1;

// above break
Point(counter)={posDamX,1.0,0, lcCoarse };
counter = counter+1;

Point(counter)={posDamX,posBreakYUp,0, lcFine };
counter = counter+1;

Point(counter)={posDamX+thicknessWall,1.0,0, lcCoarse };
counter = counter+1;

Point(counter)={posDamX+thicknessWall,posBreakYIntermUp,0, lcFine };
counter = counter+1;

Point(counter)={posDamX+thicknessWall,posBreakYUp,0, lcFine };
counter = counter+1;

// intermediate inner points
Point(counter)={posDamX+thicknessWall,0.5,0, lcFine };
counter = counter+1;

radius = 0.5*(1.0-2.0*valueIntermediate);

Point(counter)={posDamX+thicknessWall+radius,0.5,0, lcFine };
counter = counter+1;

lcounter = 1;
//+
//Line(lcounter) = {1, 5};
//lcounter = lcounter+1;
//+

//+
//Physical Curve("damwall") = {10, 9, 8, 2, 3, 4};
//+
//Physical Curve("outerwalls") = {11, 1, 7, 5};
//+
//Physical Curve("inout") = {12, 6};
//+
//+
Circle(1) = {8, 15, 13};
//+
Line(2) = {13, 12};
//+
Line(3) = {12, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {4, 7};
//+
Line(6) = {7, 8};
//+
Line(7) = {8, 9};
//+
Line(8) = {9, 6};
//+
Line(9) = {6, 5};
//+
Line(10) = {5, 1};
//+
Line(11) = {1, 2};
//+
Line(12) = {2, 10};
//+
Line(13) = {10, 11};
//+
Line(14) = {11, 14};
//+
Line(15) = {13, 14};
//+
Line(16) = {11, 6};
//+
Curve Loop(1) = {16, -8, -7, 1, 15, -14};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, 3, 4, 5, 6, 1};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {12, 13, 16, 9, 10, 11};
//+
Plane Surface(3) = {3};
//+
Physical Curve("outerwalls") = {12, 10, 3, 5};
//+
Physical Curve("inout") = {11, 4};
//+
Physical Curve("damwall") = {13, 14, 15, 2, 9, 8, 7, 6};
