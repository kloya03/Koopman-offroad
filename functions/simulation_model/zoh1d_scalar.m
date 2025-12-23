function y = zoh1d_scalar(t, t0, dt, ytab, N)
    k = (t - t0) / dt;
    tol = 1e-12;              % helps when t is 0.1 but stored as 0.10000000000001
    i = ceil(k - tol);        % makes t=t0+dt map to i=1, t=t0+2dt -> i=2, etc.

    if i < 1
        i = 1;
    elseif i > N
        i = N;
    end

    y = ytab(i);
end