function dYdt = STM_vec(t, Y, mu)

    % Compute the A matrix for the current state
    A_fi = A_CR3BP(Y(37:42), mu) * reshape(Y(1:36), [6,6]);
    A_fi_vec = reshape(A_fi, [36, 1]);

    % State derivative with the CR3BP dynamics
    dYdt = [A_fi_vec; CR3BP_mu(t, Y(37:42), mu)];
end
