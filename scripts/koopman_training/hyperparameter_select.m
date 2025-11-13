clc;
clear;

nl_sy = [400, 400, 600, 600, 600;...
    200, 300, 200, 300, 400];
nB = [200, 300, 400]; % [200, 250, 300, 350, 400];
cut_off = [2,3,4,5,6,7,8,9];

[X2, X3, X1] = ndgrid(1:numel(nB), 1:numel(cut_off), 1:size(nl_sy,2));
params = [ nl_sy(1,X1(:))', nl_sy(2,X1(:))', nB(X2(:))', cut_off(X3(:))' ];

fid = fopen('params.txt','w');
for i = 1:size(params,1)
    fprintf(fid, '%d %d %d %d\n', params(i,:));
end
fclose(fid);