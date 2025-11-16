%% 3.  find Extended subspace      %%%%%%

% Observability Matrix and system rank
function [Gam_Xi_R,rr] = find_ExObs(Xi_N,ntau)
% tau = 0.95;
[U,S,~] = svd(Xi_N,'econ');     % S is diag of singular values
sig = diag(S);                  % for PSD, sig == eigenvalues
Cu  = cumsum(sig) / sum(sig);   % NO sqrt here
rr   = max(find(round(Cu,ntau) == 1,1, 'first'),3);
ct = 0;
while isempty(rr)
    ct =ct +1;
    rr = max(find(round(Cu,ntau-ct) == 1,1, 'first'),3); % Default value if no rank is found
end
% rr   = max(find(Cu >= tau, 1, 'first'),3);
Gam_Xi_R = U(:,1:rr)*sqrt(S(1:rr,1:rr));

end



% ls = log10(max(diag(S),1e-5));
% d3 = (diff(ls));
% [~, knee] = max(d2);
% n_knee = knee + 1;