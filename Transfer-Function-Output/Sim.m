function O = Sim(ClosedLoop_z,Input,Output,k)

% This function calculates the Output of a transfer function at instant of time k using difference
% equation.
% Inputs: 
%     ClosedLoop_z  =  The transfer function in the form of
%                                         b0 + b1 z^-1 + b2 z^-2 +...+ bnb z^-nb
%                                        -----------------------------------------
%                                          1 + a1 z^-1 + a2 z^-2 +...+ ana z^-na
%     Input  =  A vector of the Inputs till the instant k
%     Output  =  A vector of the Outputs till the instant k-1
%     k  =  The instant at which Output is required
% Outputs:
%     O  =  The Output calculated at instant k

[num,den]=tfdata(ClosedLoop_z);
AA = cell2mat(den);
A = AA./AA(1);    % to make sure that the first term is 1
B = cell2mat(num);
B = B./AA(1);
delay = get(ClosedLoop_z,'iodelay');
if delay >0
B = [zeros(1,delay),B];
end
% To get delay "d ==  num of first zeros in Bz"
ind = find(B == 0);
test = isempty(ind);
if test == 1
        d = 0;
elseif ind(1) == 1
        d = 1;
elseif ind(1)~= 1 
        d = 0;
end
if length(ind)>1
for i=1:length(ind)-1
    if ind(i+1)-ind(i) == 1
       d = i+1;
    else
        break
    end
end
end

if length(A)>1
     A = A(2:end);
end
B = B(d+1:end);
na = length(A);
nb = length(B);
    %y(k)=-A(1)*y(k-1)-A(2)*y(k-2)...-A(na)*y(k-na) +B(1)*U(k)+B(2)*U(k-1)...B(nb)*U(k-nb);
    for n = 1:na
        if k-n>0
            g(n) = -A(n)*Output(k-n);
        else
            g(n) = 0;
        end
    end
    m=0:1:nb;
    for v =1:nb
        if k-d-m(v)>0
            h(v) = B(v)*Input(k-d-m(v));
        else
            h(v) = 0;
        end
    end
    O = sum(g)+sum(h);

end