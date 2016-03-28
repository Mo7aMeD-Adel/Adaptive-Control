function [a,b] = BatchLeastSquares(u,y,d,nb,na)
% This function Identify the system parameters for a known system input and
% output.
% The system Transfer Function is as in the following form:
%
%         z^-d * (bo + b1*z^-1 + b2*z^-2 + ... + b_nb*z^-nb)
% G(z) =  ---------------------------------------------------
%                1 + a1*z^-1 + ... + a_na*z^-na
%
% OUTPUTS:
% a and b are vectors of the system estimated parameters.
% INPUTS:
% u: is the system input raw vector, y is the system output raw vector.
% d: is the delay
% nb: is the number of zeros of the equired system model.
% na: is the number of poles of the equired system model.

ne = length(y);
nu = na+nb+1;

if ne == nu || nu > ne
   disp('Error: Number of Unknowns >= Number of Equations. ') 
end

for i = 1:ne
    for j = 1:nu
        if j <= na % terms of y
            if (i-j)<=0
                psi(i,j) = 0;
            else
                psi(i,j) = -y(i-j);
            end
        else       % terms of u
            if (i-d-(j-(na+1)))<=0
                psi(i,j) = 0;
            else
                psi(i,j) = u(i-d-(j-(na+1)));
            end
        end
    end
end

s = size(psi);
if s(1) == ne && s(2) == nu
    
theta = inv(psi'*psi)*psi'*y';
a = theta(1:na);
b = theta(na+1:end);

else
    disp('ERROR in epsi')
end

end