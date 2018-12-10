// Gmsh project created on Fri Dec  7 09:45:02 2018

// coarseness definition
lcCoarse = 0.03; // 0.05
lcFine = 0.005;

counter = 1;

// Geometry settings
posDamX = 0.5;
posBreakYLow = 0.35;
posBreakYUp = posBreakYLow+0.3;
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

Point(counter)={posDamX+thicknessWall,posBreakYLow,0, lcFine };
counter = counter+1;

// above break
Point(counter)={posDamX,1.0,0, lcCoarse };
counter = counter+1;

Point(counter)={posDamX,posBreakYUp,0, lcFine };
counter = counter+1;

Point(counter)={posDamX+thicknessWall,1.0,0, lcCoarse };
counter = counter+1;

Point(counter)={posDamX+thicknessWall,posBreakYUp,0, lcFine };
counter = counter+1;

lcounter = 1;
//+
Line(lcounter) = {1, 5};
lcounter = lcounter+1;
//+
Line(lcounter) = {5, 6};
lcounter = lcounter+1;
//+
Line(lcounter) = {6, 8};
lcounter = lcounter+1;
//+
Line(lcounter) = {8, 7};
lcounter = lcounter+1;
//+
Line(lcounter) = {7, 4};
lcounter = lcounter+1;
//+
Line(lcounter) = {4, 3};
lcounter = lcounter+1;
//+
Line(lcounter) = {3, 11};
lcounter = lcounter+1;
//+
Line(lcounter) = {11, 12};
lcounter = lcounter+1;
//+
Line(lcounter) = {12, 10};
lcounter = lcounter+1;
//+
Line(lcounter) = {10, 9};
lcounter = lcounter+1;
//+
Line(lcounter) = {9, 2};
lcounter = lcounter+1;
//+
Line(lcounter) = {2, 1};
lcounter = lcounter+1;
//+
Curve Loop(1) = {12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
//+
Plane Surface(1) = {1};
//+
Physical Curve("damwall") = {10, 9, 8, 2, 3, 4};
//+
Physical Curve("outerwalls") = {11, 1, 7, 5};
//+
Physical Curve("inout") = {12, 6};
