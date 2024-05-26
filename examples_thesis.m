% Examples from chapter 4 of the thesis

d = 1.0e-12;

A = triu(ones(20,20));
F20 = gallery('frank',20);
T20 = gallery('sampling',20);
J100 = gallery('jordbloc',100,0);

[N_A, D_A, R_A, A_new, sweep_A] = closest_normal(A,d);
[N_F20, D_F20, R_F20, F20_new, sweep_F20] = closest_normal(F20, d);
[N_T20, D_T20, R_T20, T20_new, sweep_T20] = closest_normal(T20,d);
[N_J100, D_J100, R_J100, J100_new, sweep_J100] = closest_normal(J100,d);

% we recreate figure 4.1.(a) in the thesis
m = 150; % number of iterations we are interested in looking at
L = zeros(1,m);
L0 = norm(diag(diag(A)),"fro");
for i=1:m
    [N, D, R, F] = closest_normal_wsweep(A,i);
    L(1,i) = norm(D,"fro");
end 

excerptL = [L0 L]; 

figure
plot(reshape(0:m, 1, (m+1)),excerptL,'-','Linewidth',2,'Color','#800000','Marker','o');

% we recreate figure 4.1.(b)
m2 = 150; % number of iterations we are interested in looking at
Diff = zeros(1,m2);
for i=1:m2
    Diff(1,i) = excerptL(1,i+1)-excerptL(1,i);
end

figure
semilogy(reshape(1:m2, 1, m2),Diff,'-','Linewidth',2,'Color','#800000','Marker','o');


% we recreate table 4.2
diffsizes = [12 20 50];
tolerances = [1.0e-4 1.0e-6 1.0e-8 1.0e-10];
sweeps = zeros(size(diffsizes,2),size(tolerances,2));
for j=1:size(diffsizes,2)
    for i=1:size(tolerances,2)
        [N, D, R, A_new, sweep_A] = closest_normal(gallery('frank',diffsizes(:,j)),tolerances(:,i));
        sweeps(j,i) = sweep_A;
    end
end


