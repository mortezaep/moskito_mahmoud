distance = 200;
depth = -6000;
ms = 20;

Point(1) = {distance/2, 0, 0, ms};
Point(2) = {distance/2, 0, -200, ms};
Point(3) = {distance/2, 0, -3500, ms};
Point(4) = {500, 250, -5000, ms};
Point(5) = {950, 600, -5900, ms};
Point(6) = {1000, 1000, depth, ms};
Point(7) = {1000, 2000, depth, ms};
Point(8) = {700, 2990, depth, ms};
Point(9) = {0, 3000, depth, ms};
Point(10) = {-700, 2990, depth, ms};
Point(11) = {-1000, 2000, depth, ms};
Point(12) = {-1000, 1000, depth, ms};
Point(13) = {-950, 600, -5900, ms};
Point(14) = {-500, 250, -5000, ms};
Point(15) = {-distance/2, 0, -3500, ms};
Point(16) = {-distance/2, 0, -200, ms};
Point(17) = {-distance/2, 0, 0, ms};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
BSpline(4) = {4, 5, 6};
Line(5) = {6, 7};
BSpline(6) = {7, 8, 9};
BSpline(7) = {9, 10, 11};
Line(8) = {11, 12};
BSpline(9) = {12, 13, 14};
Line(10) = {14, 15};
Line(11) = {15, 16};
Line(12) = {16, 17};

Physical Point("inlet") = {1};
Physical Point("outlet") = {17};

Physical Line("0_200_i") = {1};
Physical Line("200_3500_i") = {2};
Physical Line("3500_5000_i") = {3};
Physical Line("not_cased_i") = {4,5,6};
Physical Line("not_cased_p") = {7,8,9};
Physical Line("3500_5000_p") = {10};
Physical Line("200_3500_p") = {11};
Physical Line("0_200_p") = {12};