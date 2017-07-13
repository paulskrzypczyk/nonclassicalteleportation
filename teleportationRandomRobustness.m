function [TRR,Fax,G] = teleportationRandomRobustness(sigax,psix,varargin)
%TELEPORTATIONRANDOMROBUSTNESS calculates the teleportation robudstness of 
% a set of teleportation data
% This function has two required arguments:
%  sigax: a 4-D array, containing the unnormalised teleported states. The 
%  first two dimensions contain the unnormalised quantum states, while the
%  remaining two dimensions are (a,x), such that sigma(:,:,a,x) =
%  \sigma_a|psi_x.
%  psix: a 3-D array, containing the unknown states given to Alice. The
%  first two dimensions contain the unknown states, while the remaining
%  dimension is x, such that psix(:,:,x) = psi_x;
%
% TRR = teleportationRandomRobustness(sigax,psix) returns the random
% teleportation robustness of the teleportation data sigax and psix.
%
% [TRR, Fax,G] = teleportationRandomRobustness(sigax,psix) also returns the 
% teleportation witness {Fax,G} that certifies that the teleportation 
% robustness is TRR. Fax is a 4-D array, with the first two dimensions 
% containing members of the teleportation witness that multiply sig_a|psi_x
% and the last two labelling (a,x). G is a 2-D array containing the member
% of the teleportation witness that multiplies rho^B.
%
% This function has one optional argument:
%   k: (default 1)
%
% [TRR, Fax] = teleportationRandomRobustness(sigax,psix,k) calculates the
% teleportation robustness demanding that the operators MaVB have a
% k-symmetric PPT extension, instead of demanding only that it is PPT.
%
%   requires: CVX (http://cvxr.com/cvx/), QETLAB (http://www.qetlab.com)
%   authors: Paul Skrzypczyk, Daniel Cavalcanti and Ivan Supic
%   last updated: July 13, 2017

% set optional argument defaults: K=1
[k] = opt_args({1},varargin{:});

[dV,~,ma] = size(psix); % dimensions input states, and number of states
[dB,~,oa,~] = size(sigax); % dimension of Bob, number of outcomes for Alice


rhoB = squeeze(sum(sigax(:,:,:,1),3));
% reduced density operator of Bob

cvx_begin quiet

    variable MaVB(dV*dB,dV*dB,oa) hermitian semidefinite
    variable pax(oa,ma) nonnegative
    
    dual variables F{oa,ma}
    dual variable G
    
    minimise sum(sum(pax))/ma
    % 1/|x|*sum_a,x p(a|psi_x)
    
    subject to
    
    for x = 1:ma
        for a = 1:oa
            F{a,x} : sigax(:,:,a,x) + pax(a,x)*eye(dB)/dB ...
                == PartialTrace(MaVB(:,:,a)*Tensor(psix(:,:,x),eye(dB)),1,[dV,dB])
        end
    end
    % sig_a|psi_x + p(a|psi_x) id/d == tr_V[M_a^VB(psi_x otimes id)]
    
    for a = 1:oa
        SymmetricExtension(MaVB(:,:,a),k,[dV,dB],1,1) == 1;
        % M_a^VB should have a k-symmetric PPT extension (relaxation of the
        % separability requirement)
    end
    
    G : sum(MaVB,3) == Tensor(eye(dV),rhoB + sum(sum(pax))*eye(dB)/dB/ma);
    % sum_a M_a^VB == id otimes (rho^B + 1/|x| sum_a,x p(a|psi_x) id/d) 
    
cvx_end

TRR = cvx_optval;

Fax = zeros(dB,dB,oa,ma);
for a = 1:oa
    for x = 1:ma
        Fax(:,:,a,x) = -F{a,x};
    end
end

G = -G;

return