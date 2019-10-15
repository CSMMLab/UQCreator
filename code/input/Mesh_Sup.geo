lcFine = 0.03;
lcFN = 0.07;
lcNF = 0.12;
lcNormal = 0.15;
Point(1) = {0, 1, 0, lcFine};
Point(2) = {0.451145, 1.17431, 0, lcFine};
Point(3) = {1.5374, 1.58609, 0, lcFine};
Point(4) = {1.84392, 1.69785, 0, lcFine};
Point(5) = {2.09192, 1.78471, 0, lcFine};
Point(6) = {2.31702, 1.86036, 0, lcFine};
Point(7) = {2.53868, 1.93173, 0, lcFine};
Point(8) = {2.75863, 1.99947, 0, lcFine};
Point(9) = {2.98017, 2.06463, 0, lcFine};
Point(10) = {3.20553, 2.12782, 0, lcFine};
Point(11) = {3.43637, 2.18938, 0, lcFine};
Point(12) = {3.67404, 2.24954, 0, lcFine};
Point(13) = {3.91967, 2.30841, 0, lcFine};
Point(14) = {4.17432, 2.36602, 0, lcFine};
Point(15) = {4.43896, 2.42236, 0, lcFine};
Point(16) = {4.71452, 2.47738, 0, lcFine};
Point(17) = {5.00193, 2.53097, 0, lcFine};
Point(18) = {5.30215, 2.58301, 0, lcFine};
Point(19) = {5.6161, 2.63333, 0, lcFine};
Point(20) = {5.9447, 2.68172, 0, lcFine};
Point(21) = {6.28903, 2.72796, 0, lcNormal};
Point(22) = {6.65014, 2.77179, 0, lcNormal};
Point(23) = {7.02916, 2.8129, 0, lcNormal};
Point(24) = {7.42641, 2.8509, 0, lcNormal};
Point(25) = {7.84468, 2.88554, 0, lcNormal};
Point(26) = {8.28457, 2.91635, 0, lcNormal};
Point(27) = {8.74748, 2.94286, 0, lcNormal};
Point(28) = {9.23488, 2.96457, 0, lcNormal};
Point(29) = {9.74835, 2.98089, 0, lcNormal};
Point(30) = {10.2896, 2.99121, 0, lcNormal};
Point(31) = {10.8604, 2.99484, 0, lcNormal};
Point(32) = {10.8604, 0, 0, lcNormal};
Point(33) = {0, 0, 0, lcFine};
Point(34) = {0, 1, 0, lcNormal};

// lc points on center line
Point(105) = {4, -0, 0, lcFine};
Point(106) = {7, 0, 0, lcFN};
Point(107) = {8.8, -0, 0, lcNF};
Point(108) = {2, -0, -0, lcFine};

// center line
Line(81) = {33, 108};
Line(82) = {108, 105};
Line(83) = {105, 106};
Line(84) = {106, 107};
Line(85) = {107, 32};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 17};
Line(17) = {17, 18};
Line(18) = {18, 19};
Line(19) = {19, 20};
Line(20) = {20, 21};
Line(21) = {21, 22};
Line(22) = {22, 23};
Line(23) = {23, 24};
Line(24) = {24, 25};
Line(25) = {25, 26};
Line(26) = {26, 27};
Line(27) = {27, 28};
Line(28) = {28, 29};
Line(29) = {29, 30};
Line(30) = {30, 31};
Line(31) = {31, 32};
//Line(32) = {32, 33};
Line(33) = {33, 34};
Line(34) = {34, 1};
//+
Point(35) = {-1, 1, 0, lcFine};
//+
Point(36) = {-1, 0, 0, lcFine};
//+
Point(37) = {-1, 2, 0, lcNormal};
//+
Point(38) = {-3, 2, 0, lcNormal};
//+
Point(39) = {-3, 0, 0, lcNormal};
//+
Line(35) = {1, 35};
//+
Line(36) = {35, 37};
//+
Line(37) = {37, 38};
//+
Line(38) = {38, 39};
//+
Line(39) = {39, 36};
//+
Line(40) = {36, 33};
//+
Line(41) = {35, 36};
//+
//Curve Loop(36) = {38, 39, -41, 36, 37};
//+
//Plane Surface(37) = {36};
//+
//Curve Loop(37) = {41, 40, 33, 34, 35};
//+
//Plane Surface(38) = {37};

Symmetry {0, 1, 0, 0} {
  Duplicata { Curve{38}; Point{38}; Point{39}; Point{37}; Point{35}; Point{36}; Curve{33}; Curve{1}; Point{33}; Point{1}; Point{2}; Point{3}; Point{4}; Curve{4}; Point{5}; Point{6}; Point{7}; Point{8}; Curve{8}; Point{9}; Point{10}; Point{11}; Point{12}; Point{13}; Point{14}; Point{15}; Point{16}; Point{17}; Point{18}; Curve{18}; Point{19}; Point{20}; Point{21}; Point{22}; Curve{23}; Point{23}; Point{24}; Point{25}; Point{26}; Point{27}; Point{28}; Point{29}; Point{30}; Point{31}; Point{32}; Curve{2}; Curve{31}; Curve{30}; Curve{29}; Curve{28}; Curve{27}; Curve{26}; Curve{25}; Curve{24}; Curve{22}; Curve{21}; Curve{20}; Curve{19}; Curve{17}; Curve{16}; Curve{15}; Curve{14}; Curve{13}; Curve{12}; Curve{11}; Curve{10}; Curve{9}; Curve{7}; Curve{6}; Curve{5}; Curve{3}; Curve{39}; Curve{37}; Curve{36}; Curve{41}; Curve{40}; Curve{35}; Point{106}; Point{107}; Point{108}; Curve{81}; Curve{82}; Curve{83}; Curve{84}; Curve{85}; }
}

//+
Characteristic Length {119} = lcFine;


//+
Curve Loop(1) = {37, 38, 39, -41, 36};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {86, 39, -122, 121, 120};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {35, 41, 40, 33};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {81, 82, 83, 84, 85, -31, -30, -29, -28, -27, -26, -25, -24, -23, -22, -21, -20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, -33};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {81, 82, 83, 84, 85, -94, -95, -96, -97, -98, -99, -100, -101, -92, -102, -103, -104, -105, -91, -106, -107, -108, -109, -110, -111, -112, -113, -114, -90, -115, -116, -117, -89, -118, -93, -88, -87};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {124, 122, 40, 87};
//+
Plane Surface(6) = {6};

//+
Physical Curve("walls") = {30, 29, 28, 27, 26, 25, 24, 22, 23, 21, 20, 19, 17, 18, 16, 15, 13, 12, 11, 10, 8, 6, 4, 2, 1, 35, 36, 37, 120, 121, 124, 88, 93, 118, 89, 117, 116, 115, 114, 90, 113, 112, 111, 110, 109, 108, 107, 106, 91, 105, 104, 103, 92, 102, 101, 100, 99, 98, 97, 96, 95, 3, 5, 7, 9, 14, 38, 86};
//+
Physical Curve("outlet") = {31, 94};
