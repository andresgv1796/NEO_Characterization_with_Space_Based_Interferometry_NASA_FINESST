function V = sample_simulated_visibility(simulated, u, v)
%SAMPLE_SIMULATED_VISIBILITY
%
% Evaluate the Fourier transform of the simulated brightness distribution
% at nonuniform UV coordinates.
%
% Inputs:
%   simulated.I     : normalized sky brightness image
%   simulated.L_rad : l-coordinate grid [rad]
%   simulated.M_rad : m-coordinate grid [rad]
%   u, v            : UV coordinates in wavelengths
%
% Output:
%   V : complex visibility samples
%
% Convention:
%   V(u,v) = sum_p I_p exp[-2*pi*i*(u*l_p + v*m_p)]

    u = u(:);
    v = v(:);

    I = simulated.I(:);
    l = simulated.L_rad(:);
    m = simulated.M_rad(:);

    nUv = numel(u);
    V = zeros(nUv, 1);

    parfor k = 1:nUv
        phase = -2*pi*(u(k)*l + v(k)*m);
        V(k) = sum(I .* exp(1i*phase));
    end
end