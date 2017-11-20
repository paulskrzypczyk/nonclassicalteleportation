function sigax = genTeleportationData(rho,Ma,psix)
%GENTELEPORTATIONDATA outputs the teleportation data that arises from a
%specific choice of shared state, measurement, and unknown input states
% This function has three required arguments:
%  rho: a 2-D array, containing the dAdB x dAdB density operator shared by
%  Alice and Bob
%  Ma: a 3-D array, containing the oa dVdA x dVdA POVM elements of the
%  joint measurement performed by Alice, on the unknown input state, and
%  her half of the shared state with Bob. The first two dimensions contain
%  the POVM elements, whilst the last contains the outcome a.
%  psix: a 3-D array, containing the ma dV x dV density operators of the
%  unknown input states given to Alice. The first two dimensions contain
%  the states, whilst the last labels the input state x. 
%
% sigax = genTeleportationData(rho,Ma,psix) returns the unnormalised states
% prepared for Bob, after Alice applies the measurement to her half of the
% shared state, and the unknown input state. sigax is a 4-D array. The
% first two dimensions contain the dB x dB (unnormalised) density operators 
% of Bob. The third and fourth dimensions are (a,x), the outcome of Alice,
% and the label of the unknown state. 
%
%   requires: QETLAB (http://www.qetlab.com)
%   authors: Paul Skrzypczyk, Daniel Cavalcanti and Ivan Supic
%   last updated: July 13, 2017

[dV,~,ma] = size(psix); % dimension of unknown states, and number of them
[dAdB,~] = size(rho); % dimension of AB
[dVdA,~,oa] = size(Ma); % dimension of VA, and number of outcomes of Ma
dA = dVdA/dV; 
dB = dAdB/dA;

sigax = zeros(dB,dB,oa,ma);

for x = 1:ma
    for a = 1:oa
        sigax(:,:,a,x) = PartialTrace(Tensor(Ma(:,:,a),eye(dB))...
            *Tensor(psix(:,:,x),rho),[1,2],[dV,dA,dB]);
    end
end
% sig_a|psi_x = tr_VA[(M_a^VA otimes id^B)(psi_x^V otimes rho^AB)]

return