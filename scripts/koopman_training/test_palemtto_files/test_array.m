% test_array.m
disp('--- MATLAB Test Script Running ---');
ver
% Display incoming parameters
disp(['Parameter a = ', num2str(a)]);
disp(['Parameter b = ', num2str(b)]);
disp(['param_tag = ', param_tag]);

% Check SLURM IDs
job_id  = getenv('SLURM_JOB_ID');
task_id = getenv('SLURM_ARRAY_TASK_ID');

if isempty(job_id),  job_id  = 'local'; end
if isempty(task_id), task_id = '0';     end

disp(['SLURM job_id: ', job_id]);
disp(['SLURM task_id: ', task_id]);

% Save workspace file
ws_name = sprintf('abc/ws_test_%s_job%s_task%s.mat', param_tag, job_id, task_id);

param_sum = a+b
save(ws_name);


disp(['Saved workspace as: ', ws_name]);
disp('--- Done ---');
