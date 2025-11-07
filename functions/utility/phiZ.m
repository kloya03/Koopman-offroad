%% 8. PhiZ
function [phi_Z_N,D_PhiZ,Ex_obs] = phiZ(Cr,Kr,nB)      %%%%%%
rr = size(Kr,1);
ny = size(Cr,1);
Ex_obs = zeros(ny,rr,nB+1);
phi_Z = zeros(ny*nB,rr);
for ii = 1:nB
    Ex_obs(:,:,ii+1) =  Cr*(Kr^(ii-1));
    phi_Z(ny*(ii-1)+1:ny*ii,:) = Cr*(Kr^(ii-1));
end
D_PhiZ = std(phi_Z, 0, 1);
phi_Z_N = (phi_Z)./D_PhiZ;% (phi_Z - mean(phi_Z, 1))./std(phi_Z, 0, 1); % normalize
end
