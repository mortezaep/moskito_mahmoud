ms = 10;

Point(1) = {0, 0, 0, ms};
Point(2) = {0, 0, -500, ms};

Line(1) = {1, 2};

Physical Point("inlet") = {1};
Physical Point("outlet") = {2};

Physical Line("well") = {1};

