distance = 200;
depth = -6000;
ms = 20;

Point(1) = {distance/2, 0, 0, ms};
Point(2) = {distance/2, 0, -3500, ms};
Point(3) = {500, 250, -5000, ms};
Point(4) = {950, 600, -5900, ms};
Point(5) = {1000, 1000, depth, ms};
Point(6) = {1000, 2000, depth, ms};
Point(7) = {700, 2990, depth, ms};
Point(8) = {0, 3000, depth, ms};
Point(9) = {-700, 2990, depth, ms};
Point(10) = {-1000, 2000, depth, ms};
Point(11) = {-1000, 1000, depth, ms};
Point(12) = {-950, 600, -5900, ms};
Point(13) = {-500, 250, -5000, ms};
Point(14) = {-distance/2, 0, -3500, ms};
Point(15) = {-distance/2, 0, 0, ms};

Line(1) = {1, 2};
Line(2) = {2, 3};
BSpline(3) = {3, 4, 5};
Line(4) = {5, 6};
BSpline(5) = {6, 7, 8};
BSpline(6) = {8, 9, 10};
Line(7) = {10, 11};
BSpline(8) = {11, 12, 13};
Line(9) = {13, 14};
Line(10) = {14, 15};

Physical Point("inlet") = {1};
Physical Point("outlet") = {15};

Physical Line("0_2000_i") = {1};
Physical Line("2000_3500_i") = {2};
Physical Line("3500_4500_i") = {3};
Physical Line("horizontal_i") = {4};
Physical Line("curved_i") = {5};
Physical Line("curved_p") = {6};
Physical Line("horizontal_p") = {7};
Physical Line("3500_4500_p") = {8};
Physical Line("2000_3500_p") = {9};
Physical Line("0_2000_p") = {10};