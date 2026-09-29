function [EGMOS, delta] = gmos_accuracy_test_paper(z, N, mu, ode_opts)
%GMOS_ACCURACY_TEST_PAPER Paper-style GMOS accuracy test (Eq. 25).
%
% Inputs
%   z       : packed solution [X0(:); T; rho]
%   N       : number of curve points
%   mu      : CR3BP mass parameter
%   ode_opts: odeset options
%
% Outputs
%   EGMOS : scalar paper accuracy metric
%   delta : 1xN vector of per-point errors

    if nargin < 4 || isempty(ode_opts)
        ode_opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    end

    [X0, T, rho] = gmos_unpack_z(z, N);

    % Midpoint rotation alpha = pi/N (paper)
    alpha = pi / N;
    Xmid = gmos_shift_curve_fft(X0, -alpha);   % X(theta + alpha) = X(theta - (-alpha))

    % Propagate midpoints one strobe time
    Xmid_1 = gmos_flowmap_curve_CR3BP(Xmid, T, mu, ode_opts);

    % Rotate back by rho
    Xmid_1_shift = gmos_shift_curve_fft(Xmid_1, rho);

    % delta_i = || R_-rho(Psi_T(v(theta_i+alpha))) - v(theta_i+alpha) ||
    D = Xmid_1_shift - Xmid;
    delta = sqrt(sum(D.^2, 1));

    % Eq. (25)
    EGMOS = sqrt(mean(delta.^2));
end