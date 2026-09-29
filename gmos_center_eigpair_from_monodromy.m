function cen = gmos_center_eigpair_from_monodromy(tspan, xPO0, mu, mode)
%GMOS_CENTER_EIGPAIR_FROM_MONODROMY Select a GMOS center eigenpair.
%
% Thin production wrapper:
%
%   1. compute the one-period CR3BP monodromy eigensystem;
%   2. delegate ALL center-mode selection to
%      gmos_center_eigpair_from_eigs.
%
% Keeping the actual selection rule in gmos_center_eigpair_from_eigs avoids
% duplicating trivial-pair removal, unit-circle tolerances, and planar /
% vertical mode-selection logic in analysis scripts.

    if nargin < 4 || isempty(mode)
        mode = 'planar';
    end

    M_eigs = monodromy_CR3BP(tspan,xPO0,mu);

    V = M_eigs(:,1:6);
    D = M_eigs(:,7:12);

    cen = gmos_center_eigpair_from_eigs(V,D,mode);

end
