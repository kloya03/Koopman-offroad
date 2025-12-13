function [U,X_true,cost,ff] = MPC_koopman(A,B,C,x0,yref,T,str,CV)
%%%%%%%%%%% MODEL PREDICTIVE CONTROL %%%%%%%%%%%
load(str,"Data")
f_u = Data.f_u; f_ud = Data.f_ud; const = Data.const;
tsim = T(:,1); tpred = T(:,2); tcont =T(:,3);
dt = T(:,4); tc=0; nsim = 0;
xc = x0; eps = 1e-10;
fal = []; f_act = [];
U=[]; Xtru = xc; uminmax=3;
r = size(A,1);        % number of states
nu = size(B,2);        % number of inputs
ny = size(C,1);     % number of outputs
Q = CV(1:ny,1); QN = CV(1:ny,2); R = CV(1:nu,3);
%%
load(str,"MDL_fit")
Z_GP =zeros(r,2);
for i=1:r
    [Zi_mean,Zi_sd] = predict(MDL_fit(i).gprMdl,x0.');
    Z_GP(i,:) = [Zi_mean Zi_sd.^2];
end
Z0 = Z_GP(:,1);
x0_act=C*Z0;
abs(x0_act - x0)

%% MPC Begins
tic
while tc < tsim -eps
    if tsim-tc <=tcont
        Np = round((tsim-tc)/dt);
        Nc = Np;
    elseif (tcont < tsim-tc)&& (tsim-tc <= tpred)
        Np = round((tsim-tc)/dt);
        Nc = round(tcont/dt);
    else
        Np = round(tpred/dt);
        Nc = round(tcont/dt);
    end
    ri = R.*eye(nu); rNp = R.*eye(nu);
    qi = Q.*eye(ny); qNp = QN.*eye(ny);
    Yr = [yref(:,round(tc/dt)+(1:Np))];
    lb = -uminmax*ones(nu*(Np),1);
    ub = uminmax*ones(nu*(Np),1);

    [H,f,Phi,Gamma] = MPC_matrices_noD(Z0,A,B,C,Yr,ri,rNp,qi,qNp,Np);
    Hbar = (H+H.')/2;
    options = optimoptions('quadprog','Display','iter-detailed' );
    [Uqp,fval,exitflag,output] = quadprog(Hbar,f,[],[],[],[],...
        lb,ub,[],options);

    % Trajectory Extraction
    fal = [fal fval];
    u_GP = [Uqp(1:nu*Nc,1)];
    U = [U;Uqp(1:nu*Nc,1)];

    nsim = nsim+Nc;
    Xtr = xc;
    for jj=1:Nc
        Xnext = f_ud(nsim+(jj*dt),Xtr(:,jj),const,u_GP(jj,1));
        Xtr = [Xtr, Xnext];
        Xtru = [Xtru, Xnext];
    end
    Zpred = reshape(Phi*Z0 + Gamma*Uqp,r,Np);
    ypred = C*Zpred;
    f_act = [f_act sum(diag((ypred(:,1:end-1)-...
        Yr(:,1:end-1)).'*qi*(ypred(:,1:end-1)-Yr(:,1:end-1))))...
        + sum(diag(Uqp.'*ri*Uqp)) + (ypred(:,end)-...
        Yr(:,end)).'*qNp*(ypred(:,end)-Yr(:,end))];

    xc = Xtr(:,end);

    for i=1:r
        [Zi_mean] = predict(MDL_fit(i).gprMdl,xc.');
        Z0(i,1) = Zi_mean;
    end

    tc = round(tc + Nc*dt,3)
end

ff= [f_act.', fal.'];
X_true = x0;
for jj=1:nsim
    Xnext = f_ud(jj*dt,X_true(:,jj),const,U(jj,1));
    X_true = [X_true, Xnext];
end
cost = dt*CostFcn(X_true.',U,CV,yref);

end
