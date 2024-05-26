% Examples from Ruhe's paper, recalculated using our implementation

d = 10*eps;

A = [0.7616+1.2296i -1.4740-0.4577i; -1.6290-2.6378i 0.1885-0.8575i];
B7 = rand(10,10);
J7 = gallery('jordbloc',7,0);
F12 = gallery('frank',12);

[N_A, D_A, R_A, A_new, sweep_A] = closest_normal(A, d);
[N_B7, D_B7, R_B7, B7_new, sweep_B7] = closest_normal(B7, d);
[N_J7, D_J7, R_J7, J7_new, sweep_J7] = closest_normal(J7, d);
[N_F12, D_F12, R_F12, F12_new, sweep_F12] = closest_normal(F12, d);

Y_A = obesity_Y(A_new);
Y_B7 = obesity_Y(B7_new);
Y_J7 = obesity_Y(J7_new);
Y_F12 = obesity_Y(F12_new);

spread_Y_A = spread(Y_A);
spread_Y_B7 = spread(Y_B7);
spread_Y_J7 = spread(Y_J7);
spread_Y_F12 = spread(Y_F12);

Examples_Ruhe = ["A"; "B7"; "J7"; "F12"];
Spreads_Ruhe = [spread_Y_A; spread_Y_B7; spread_Y_J7; spread_Y_F12];
results_Ruhe_examples = table(Examples_Ruhe, Spreads_Ruhe) % a small table comparing the spread values calculated for A, B7, J7, and F12