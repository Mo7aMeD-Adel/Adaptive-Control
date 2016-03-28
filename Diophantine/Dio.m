function [R,S] = Dio(A,B,Alpha,d)
% Diophantine Equation Solver

% Note: Polynomials are in Z^-1 form.

% A : [1,  a_1,  a_2,  a_3,  ..., a_na] 
% B : [b_0,  b_1, b _2,  b_3,  ..., a_nb] 
% R : [1,  r_1,  r_2,  r_3,  ..., r_nr] 
% S : [s_0,  s_1, s _2,  s_3,  ..., s_ns]
% Alpha : [1, Alpha_1, Alpha_2, ..., Alpha_n_alpha] 
% d: System delay
%
%                       z^-d * B(z^-1) S(z^-1)
% G_cl(z^-1) =  -----------------------------------------
%                A(z^-1) R(z^-1) + z^-d  B(z^-1) S(z^-1)

%% Calculation of Orders

na = length(A)-1;
nb = length(B)-1;
nr = nb+d-1;
ns = na-1;
n_alpha = na+nb+d-1;

%% Variables Initialization
X = zeros(n_alpha,n_alpha);

%% Solution
for m = 1:n_alpha
    if m<=nr
        X(m:na+m,m)= A;
    elseif d == 1
        X(m-nr:m-nr-1+length(B),m) = B;
    else
        Bb = [zeros(1,d-1),B];
        X(m-nr:m-nr-1+length(Bb),m) = Bb;
    end
end
F = [A(2:end),zeros(1,n_alpha-na)]';
Beta = Alpha(2:end)-F;
V = inv(X)*Beta;
R = V(1:nr);
S = V(nr+1:end);
end