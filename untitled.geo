z = -5000;
ms = 10;
x = 5000;


Point(1) = {0, 0, 0, ms};
Point(2) = {0, 0, -430, ms};
Point(3) = {0, 0, -2150, ms};
Point(4) = {0, 0, z, ms};
Point(5) = {x, 0, z, ms};
Point(6) = {x, 0, -2150, ms};
Point(7) = {x, 0, -430, ms};
Point(8) = {x, 0, 0, ms};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};


Physical Point("inlet") = {1};
Physical Point("outlet") = {8};

Physical Line("l0_430") = {1};
Physical Line("l430_2150") = {2};
Physical Line("l2150_5000") = {3};
Physical Line("horizontal") = {4};
Physical Line("r2150_5000") = {5};
Physical Line("r430_2150") = {6};
Physical Line("r0_430") = {7};

