function y = foh1d_scalar(t, t0, dt, ytab, N)
    % clamp index to [1, N-1]
    i = floor((t - t0)/dt) + 1;
    if i < 1, i = 1; elseif i >= N, i = N-1; end
    a = (t - (t0 + (i-1)*dt))/dt;       % frac in [0,1)
    y = ytab(i) + a*(ytab(i+1) - ytab(i));
end