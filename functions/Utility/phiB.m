%% 7. PhiB

function [Phi_B,D_phiB,BZ] = phiB(Ex_obs,U_tr,Y_tr,nB,Phi_Z)        %%%%%%
rr = size(Phi_Z,2);
nu = size(U_tr,1)./nB;
ny = size(Y_tr,1)./nB;

ntr = size(U_tr,2);
for i=1:ntr
    i
    UU = U_tr(:,i);
    YY = Y_tr(:,i);
    phi_B = zeros(ny*nB,nu*rr);
    for ii = 1:nB
        PHIB = zeros(ny,nu*rr);
        for kk=1:ii
            phib = kron(UU((kk-1)*nu+1:nu*kk).',Ex_obs(:,:,ii+1-kk));
            PHIB = PHIB + phib;
        end
        phi_B(ny*(ii-1)+1:ny*ii,:) =PHIB;
    end

    % Phi_B(:,:,i) = phi_B;
    D_phiB(:,:,i) = std(phi_B,0,1);
    Phi_B(:,:,i) = (phi_B)./D_phiB(:,:,i);    % normalize
    % (phi_B-mean(phi_B, 1))./std(phi_B,0,1);    % normalize
    if nargout > 2
        phi = [phi_B Phi_Z];
        BZ(:,i) = pinv((phi.')*phi)*(phi.')*YY;
    end
end
clc;

end