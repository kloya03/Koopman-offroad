function [f,g] = BZ_cost_func(x,Phi_B,Phi_Z,Yhr,nB)
    ntr = size(Phi_B,3);
    rr = size(Phi_Z,2);
    nu = size(Phi_B,2)./rr;
    ny = size(Phi_B,1)./nB;
    nn = 1/(2*nB*ntr);
    
    % Initial Guess
    Zlift0 = reshape(x(nu*rr+1:end,:),rr,ntr);  % Zlift0 = zeros(rr,ntr);  % parameter 2
    vecB = x(1:nu*rr,:);   % Br = zeros(rr,1);  % parameter 1
    
    prediction  = reshape(pagemtimes(Phi_B, vecB), ny*nB,ntr) + Phi_Z*Zlift0;   % [ny*nB x 1 x ntr] --> [ny*nB x ntr]
    residual    = Yhr - prediction;
    f = nn * sum(residual(:).^2);
    % f = f + 0.5*lambda_B*(vecB.'*vecB) + 0.5*lambda_Z*sum(Z0(:).^2);
    
    if nargout > 1
        gZ      = -2*nn*(Phi_Z.'*residual);        % [rr x ntr]
        gB_s = pagemtimes(pagetranspose(Phi_B),reshape(residual,ny*nB,1,ntr));
        gB = (-2*nn)*sum(gB_s,3);  % []
        % gB = gB + lambda_B*vecB;   gZ = gZ + lambda_Z*Z0;
    
        g = [gB; gZ(:)];
    end
end