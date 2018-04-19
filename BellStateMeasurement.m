function Ma = BellStateMeasurement(d)
%BELLSTATEMEASUREMENT generates the generalised (full) Bell State Measurement
% in dimension d
%  This function has one required arguments:
%  d: the dimension of the Hilbert space
% Ma = BellStateMeasurement(d) returns the generalised Bell State Measurement
% in dimension d as a 3-D array, of dimension d x d x d^2, where the POVM 
% elements are stored in the first two dimensions, and the third dimension
% stores the label for the individual elements.
%
%   Requires: QETLAB (http://www.qetlab.com)
%   authors: Paul Skrzypczyk, Daniel Cavalcanti and Ivan Supic
%   last updated: 19 April 2018

X = zeros(d,d); % generalised X
Z = zeros(d,d); % generalised Z

w = exp(2*pi*1i/d);

for i = 0:d-1
    for j = 0:d-1
        X(i+1,j+1) = (mod(j-i,d) == 1);
        Z(i+1,j+1) = w^i*(i == j);
    end
end

Ma = zeros(d^2,d^2,d^2);
cnt = 1;
phi = MaxEntangled(d);

for i = 0:d-1
    for j = 0:d-1
        psi = Tensor(eye(d),X^i*Z^j)*phi;
        Ma(:,:,cnt) =  psi*psi';
        cnt = cnt+1;
    end
end
