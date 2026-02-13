%% 7. PhiB

function [Phi_B,D_phiB,BZ] = phiB_parallel(Ex_obs,U_tr,Y_tr,nB,Phi_Z)        %%%%%%
rr = size(Phi_Z,2);
nu = size(U_tr,1)./nB;
ny = size(Y_tr,1)./nB;
Econst = parallel.pool.Constant(Ex_obs);
PZconst = parallel.pool.Constant(Phi_Z);
want_BZ =  nargout > 2;

ntr = size(U_tr,2);
parfor i=1:ntr
    UU = U_tr(:,i);
    YY = Y_tr(:,i);
    phi_B = zeros(ny*nB,nu*rr);
    Ex_local = Econst.Value;   % [ny x rr x nB]
    for ii = 1:nB
        PHIB = zeros(ny,nu*rr);
        for kk=1:ii
            phib = kron(UU((kk-1)*nu+1:nu*kk).',Ex_local(:,:,ii+1-kk));
            PHIB = PHIB + phib;
        end
        phi_B(ny*(ii-1)+1:ny*ii,:) =PHIB;
    end

    % Phi_B(:,:,i) = phi_B;
    % D_phiB(:,:,i) = std(phi_B,0,1);
    Phi_B(:,:,i) = (phi_B);    % normalize
    % (phi_B-mean(phi_B, 1))./std(phi_B,0,1);    % normalize
    if want_BZ
        Phi_Z_local = PZconst.Value;
        phi = [phi_B Phi_Z_local];
        BZ(:,i) = pinv((phi.')*phi)*(phi.')*YY;
    end
end

D_phiB = std(reshape(permute(Phi_B, [1 3 2]), [ny*nB*ntr, nu*rr]));
Phi_B = Phi_B./(reshape(D_phiB,1,nu*rr,1));

end