ms = 20;
d = 50;

Point(1) = {500, 0, 5000-d, ms};
Point(2) = {500, 0, 4950-d, ms};
Point(3) = {500, 0, 3000-d, ms};
Point(4) = {500, 0, 200-d, ms};
Point(5) = {500, 0, 0-d, ms};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};

Physical Point("inlet") = {5};
Physical Point("outlet") = {1};

Physical Line("5000_4950") = {1};
Physical Line("4950_3000") = {2};
Physical Line("3000_200") = {3};
Physical Line("200_0") = {4};

