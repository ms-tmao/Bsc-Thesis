% Some Jordan Matrix Examples

d = 1.0e-12;

ns1 = [2 3 5 4 2 1 3];
diagvals1 = [10 9 0 1i -1 3 5];
ns2 = [4 7 2 3 5];
diagvals2 = [7 1+1i 3 45-1i 17];

J1 = create_Jordan_matrix(ns1,diagvals1); % a Jordan matrix with diagonal 7 blocs of sizes as in ns1 and values as in diagvals1
J2 = create_Jordan_matrix(ns2,; % a Jordan matrix with diagonal 5 blocs of sizes as in ns2 and values as in diagvals2
J3 = gallery('jordbloc',30,17); % a single Jordan bloc

[N_J1, D_J1, R_J1, J_new1, sweepJ1] = closest_normal(J1, d);
[N_J2, D_J2, R_J2, J_new2, sweepJ2] = closest_normal(J2, d);
[N_J3, D_J3, R_J3, J_new3, sweepJ3] = closest_normal(J3, d);

Y_J1 = obesity_Y(J_new1);
Y_J2 = obesity_Y(J_new2);
Y_J3 = obesity_Y(J_new3);

spread_Y_J1 = spread(Y_J1);
spread_Y_J2 = spread(Y_J2);
spread_Y_J3 = spread(Y_J3);

s1 = check_G_tang(J_new1,Y_J1);
s2 = check_G_tang(J_new2,Y_J2);
s3 = check_G_tang(J_new3,Y_J3);

