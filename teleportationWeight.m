function TW = teleportationWeight(sigax,omegax,varargin)
%TELEPORTATIONWEIGHT calculates the teleportation weight of a set of 
%teleportation data
%  This function has two required arguments:
%  sigax: a 4-D array, containing the unnormalised teleported states. The 
%  first two dimensions contain the unnormalised quantum states, while the
%  remaining two dimensions are (a,x), such that sigma(:,:,a,x) =
%  \sigma_a|psi_x.
%  omegax: a 3-D array, containing the unknown states given to Alice. The
%  first two dimensions contain the unknown states, while the remaining
%  dimension is x, such that psix(:,:,x) = psi_x;
%
% TW = teleportationWeight(sigax,omegax) returns the teleportation weight
% of the teleportation data sigax and omegax.
%
% This function has one optional argument:
%   k: (default 1)
%
% TW = teleportationWeight(sigax,omegax,k) calculates the
% teleportation weight demanding that the operators MaVB have a k-symmetric PPT 
% extension, as a relaxation of separability. The default case k=1 only imposes 
% PPT as the relaxation of separability.
%
%   requires: CVX (http://cvxr.com/cvx/), QETLAB (http://www.qetlab.com)
%   authors: Paul Skrzypczyk, Daniel Cavalcanti and Ivan Supic
%   last updated: 19 April 2018

% set optional argument defaults: k=1
[k] = opt_args({1},varargin{:});


[dB,~,oa,ma] = size(sigax); % dim. of B, no. of outcomes of Ma, number of 
                            % input states for Alice
[dV,~,~] = size(omegax); % dimension of input states

X = zeros(dV^2,dV^2);
for i = 1:dV^2
    X(i,:) = reshape(omegax(:,:,i),1,[]);
end
% X will be used to tell if the set of input states is tomographically
% complete or not (needed to improve numerical stability of code)

cvx_begin sdp quiet

    variable Ma(dV*dB,dV*dB,oa) hermitian
    variable Na(dV*dB,dV*dB,oa) hermitian
    variable rhoMB(dB,dB) hermitian semidefinite
    variable rhoNB(dB,dB) hermitian semidefinite
    
    minimise trace(rhoMB)
    
    subject to
    
    for x = 1:ma
        for a = 1:oa
            sigax(:,:,a,x) == PartialTrace((Ma(:,:,a)+Na(:,:,a))*Tensor(omegax(:,:,x),eye(dB)),1,[dV,dB])
        end
    end
    % sig_a|omega_x = tr_V[(Ma + Na)(omega_x otimes Id)]
    
    sum(Ma,3) == Tensor(eye(dV),rhoMB)
    % no-signalling
    
    if rank(X) ~= dV^2
        sum(Na,3) == Tensor(eye(dV),rhoNB)
        % if the set of states is tomographically complete, this constraint is not necessary
        % and it causes problems for numerical stability to impose this redundant constaint
        % if the set of states is not tomographically complete, this is an independent constaint
    end
    
    for a = 1:oa
        SymmetricExtension(Na(:,:,a),k,[dV,dB],1,1) == 1;
        % separability of Na (approximately)
        PartialTranspose(Ma(:,:,a)) == hermitian_semidefinite(dV*dB)
        % Ma needs to be quantum realisible (i.e. PPT)
    end
    
cvx_end

TW = cvx_optval;