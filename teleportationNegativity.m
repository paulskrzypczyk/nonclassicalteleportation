function TN = teleportationNegativity(sigax,omegax)
%TELEPORTATIONNEGATIVITY provides a lower bound on the negativity of the
%quantum state used in the teleportation experiment
%  This function has two required arguments:
%  sigax: a 4-D array, containing the unnormalised teleported states. The 
%  first two dimensions contain the unnormalised quantum states, while the
%  remaining two dimensions are (a,x), such that sigma(:,:,a,x) =
%  \sigma_a|psi_x.
%  omegax: a 3-D array, containing the unknown states given to Alice. The
%  first two dimensions contain the unknown states, while the remaining
%  dimension is x, such that psix(:,:,x) = psi_x;
%
% TN = teleportationNegativity(sigax,omegax) returns the teleportation 
% negativity of the teleportation data sigax and omegax.
%
%   requires: CVX (http://cvxr.com/cvx/), QETLAB (http://www.qetlab.com)
%   authors: Paul Skrzypczyk, Daniel Cavalcanti and Ivan Supic
%   last updated: 19 April 2018

[dB,~,oa,ma] = size(sigax); % dim. of B, no. of outcomes of Ma, number of 
                            % input states for Alice
[dV,~,~] = size(omegax); % dimension of input states

cvx_begin sdp quiet

    variable Map(dV*dB,dV*dB,oa) hermitian semidefinite% M_a^VB
    variable Mam(dV*dB,dV*dB,oa) hermitian semidefinite % N_a^VB
    variable rhoBp(dB,dB) hermitian semidefinite
    variable rhoBm(dB,dB) hermitian semidefinite
    
    minimise trace(rhoBm) 
    
    subject to
    
    for x = 1:ma
        for a = 1:oa
            sigax(:,:,a,x) == PartialTrace((Map(:,:,a)-Mam(:,:,a))*...
                Tensor(omegax(:,:,x),eye(dB)),1,[dV,dB])
        end
    end
    % sig_a|x = tr_V[(Ma+ - Ma-)(omega_x otimes Id)]
    
    sum(Map,3) == Tensor(eye(dV),rhoBp)
    sum(Mam,3) == Tensor(eye(dV),rhoBm)
    % sum_a Ma+/- = Id otimes rhoB (no-signalling)
    
cvx_end

TN = cvx_optval;

