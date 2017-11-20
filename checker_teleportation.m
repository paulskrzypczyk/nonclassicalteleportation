da = 3;
db = da;

Phi = MaxEntangled(da);

w = exp(2*pi*1i/da);
Z = zeros(da,da);
X = zeros(da,da);

for i = 0:da-1
    Z(i+1,i+1) = w^i;
    X(i+1,1+mod(i+1,da)) = 1;
end


% partial BSM

oa = 2;

MaVA = zeros(da*db,da*db,2);
MaVA(:,:,1) = Phi*Phi';
MaVA(:,:,2) = eye(da*db) - MaVA(:,:,1);

% full BSM
%{

oa = da^2;

MaVA = zeros(da*db,da*db,da^2);
cnt = 1;
for i = 0:da-1
    for j = 0:da-1
        MaVA(:,:,cnt) = Tensor(Z^i*X^j,eye(da))*Phi*Phi'*Tensor(Z^i*X^j,eye(da))';
        cnt = cnt + 1;
    end
end
%}        
rho = RandomDensityMatrix(da*db);

ER = entanglementRandomRobustness(rho,2)

omegax = zeros(da,da,da^2);
for i = 1:da^2
    U = RandomUnitary(da);
    omegax(:,:,i) = U*blkdiag(1,zeros(da-1))*U';
end

sigax = genTeleportationData(rho,MaVA,omegax);

[TRR,~,~,MaVB] = teleportationRandomRobustness(sigax,omegax)

oa*ER/da^2

sumall(abs(MaVB(:,:,1) - 1/da*PartialTranspose(rho + da^2*TRR/oa*eye(da^2)/da^2,1)))

temp = zeros(da^2);
for i = 0:da-1
    for j = 0:da-1
        temp = temp + 1/da^2*Tensor(Z^i*X^j,eye(da))*PartialTranspose(rho,1)*Tensor(Z^i*X^j,eye(da))';
    end
end

Tensor(eye(da)/da,PartialTrace(rho,1)) - temp
