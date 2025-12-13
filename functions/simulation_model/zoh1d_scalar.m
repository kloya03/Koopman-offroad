function y = zoh1d_scalar(t, t0, dt, ytab, N)
    i = floor((t - t0)/dt) + 1;
    if i < 1, i = 1; elseif i > N, i = N; end
    y = ytab(i);
end