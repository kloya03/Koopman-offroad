%% 1. create Hankel Matrix    %%%%%%
function [H1,H2,H3] = createHankelMatrix(data,nl,sy,mean_std)
% H1: past hankel matrix
% H2: future hankel matrix
% H3: Hankel matrix with no repetition
% data: {1, No. of trajectory}(numPoints, numFeatures)
% mean: [mean of fetures, 1;
%        std of same features];
% nl: row size of hankel
% sy: divide into past=sy and future=nl-sy hankel
% Construct Hankel matrix using vectorized indexing
H = [];
H3 = [];
if iscell(data)
    for i=1:size(data,2)
        D = data{i};
        D_norm = ( D - mean_std(1,:) )./mean_std(2,:);
        [numPoints, numFeatures] = size(D);
        nc = numPoints-nl+1;
        nr = numFeatures*nl;
        idx = (0:nl-1)'; % Create index offsets for time steps
        HH = reshape(D_norm(idx + (1:nc), :).', nr, nc);
        H = [H, HH];
        if nargout > 2
            H3 = [H3, HH(:,1:nl:end)];
        end
    end
else
    [numPoints, numFeatures] = size(data);
    nc = numPoints-nl+1;
    nr = numFeatures*nl;
    idx = (0:nl-1)'; % Create index offsets for time steps
    data_norm = ( data - mean_std(1,:) )./mean_std(2,:);
    H = reshape(data_norm(idx + (1:nc), :).', nr, nc);
end
H1 = H(1:numFeatures*sy,:);
H2 = H(numFeatures*sy+1:end,:);
end