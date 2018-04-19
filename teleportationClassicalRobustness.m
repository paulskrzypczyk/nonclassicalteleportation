function TR = teleportationClassicalRobustness(sigax,omegax,varargin)
%TELEPORTATIONCLASSICALROBUSTNESS calculates the classical teleportation
%robustness of a set of teleportation data
%  This function has two required arguments:
%  sigax: a 4-D array, containing the unnormalised teleported states. The 
%  first two dimensions contain the unnormalised quantum states, while the
%  remaining two dimensions are (a,x), such that sigma(:,:,a,x) =
%  \sigma_a|psi_x.
%  omegax: a 3-D array, containing the unknown states given to Alice. The
%  first two dimensions contain the unknown states, while the remaining
%  dimension is x, such that psix(:,:,x) = psi_x;
%
% TR = teleportationClassicalRobustness(sigax,omegax) returns the classical
% teleportation robustness of the teleportation data sigax and omegax.
%
% This function has one optional argument:
%   k: (default 1)
%
% TR = teleportationClassicalRobustness(sigax,omegax,k) calculates the
% classical teleportation robustness demanding that the operators MaVB 
% and NaVB have a k-symmetric PPT extension, as a relaxation of separability. 
% The default case k=1 only imposes PPT as the relaxation of separability.
%
%   requires: CVX (http://cvxr.com/cvx/), QETLAB (http://www.qetlab.com)
%   authors: Paul Skrzypczyk, Daniel Cavalcanti and Ivan Supic
%   last updated: 19 April 2018

% set optional argument defaults: k=1
[k] = opt_args({1},varargin{:});


[dB,~,oa,ma] = size(sigax); % dim. of B, no. of outcomes of Ma, number of 
                            % input states for Alice
[dV,~,~] = size(omegax); % dimension of input states


cvx_begin sdp quiet

    variable Ma(dV*dB,dV*dB,oa) hermitian
    variable Na(dV*dB,dV*dB,oa) hermitian
    variable rhoB(dB,dB) hermitian semidefinite
    
    minimise trace(rhoB)
    
    subject to
    
    for x = 1:ma
        for a = 1:oa
            sigax(:,:,a,x) == PartialTrace((-Ma(:,:,a)+Na(:,:,a))*Tensor(omegax(:,:,x),eye(dB)),1,[dV,dB])
        end
    end
    % sig_a|omega_x + tr_V[Ma(omega_x otimes Id)] = tr_V[Na(omega_x otimes Id)]
    
    sum(Ma,3) == Tensor(eye(dV),rhoB)
    sum(Na,3) == Tensor(eye(dV),sum(sigax(:,:,:,1),3) + rhoB)
    % (no-signalling + consistency)
    
    for a = 1:oa
        SymmetricExtension(Na(:,:,a),k,[dV,dB],1,1) == 1;
        SymmetricExtension(Ma(:,:,a),k,[dV,dB],1,1) == 1;
    end
    % approximate separability
    
cvx_end
    

TR = cvx_optval;