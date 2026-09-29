function targetCases = simulated_target_case_library()
%SIMULATED_TARGET_CASE_LIBRARY Parameter presets for simulated target images.
%
% This function stores editable example parameter sets for
% make_simulated_target_image.m.
%
% Important conventions:
%   l_mas, m_mas  : angular offsets on the target sky plane [mas]
%   pa_deg        : position angle measured counterclockwise from +l to +m
%   flux          : relative brightness weight, not physical calibrated flux
%
% The simulated images are normalized inside make_simulated_target_image so
% that sum(I(:)) = 1. Therefore, absolute flux values only matter through
% their ratios. For example, flux2 = 0.7 means source 2 is 70 percent as
% bright as source 1 before normalization.
%
% Current-scale note:
%   These presets are tuned for the short-baseline test case with
%   rho_max ~ 1 Glambda and fov_masObs ~ 40 mas. The nominal angular
%   scale is ~0.2 mas, and for nPixObs = 1025 the pixel size is
%   ~0.078 mas.

    targetCases = struct();


    %% Delta source
    % Exact one-pixel numerical delta source for sanity checks.
    %
    % Expected behavior:
    %   centered delta     : |V| = 1, phase = 0, dirty image = dirty beam
    %   off-center delta   : |V| = 1, phase ramp, shifted dirty image

    p = struct();
    p.l0_mas = 0.0;
    p.m0_mas = 0.0;
    p.flux   = 1.0;

    targetCases.delta = p;


    %% Point source
    % Compact finite Gaussian source.
    %
    % This is not a mathematical point. It is a finite-width Gaussian used
    % to test partial resolution. With rho_max ~ 1 Glambda, sigma_mas near
    % 0.08 mas gives a compact source that is still represented on the grid.

    p = struct();
    p.l0_mas = 0.0;
    p.m0_mas = 0.0;
    p.sigma_mas = 0.08;
    p.flux = 1.0;

    targetCases.point = p;


    %% Binary source
    % Two compact sources separated by sep_mas.
    %
    % sep_mas = 0.8 mas is several times the nominal ~0.2 mas resolution
    % scale, so the binary should create clear nontrivial visibility
    % structure without being absurdly small on the image grid.

    p = struct();
    p.center_l_mas = 0.0;
    p.center_m_mas = 0.0;
    p.sep_mas = 0.80;
    p.pa_deg = 20.0;
    p.sigma1_mas = 0.08;
    p.sigma2_mas = 0.08;
    p.flux1 = 1.0;
    p.flux2 = 0.7;

    targetCases.binary = p;


    %% Star + planet
    % Bright compact central source plus a fainter off-axis companion.
    %
    % This is intentionally high contrast only by real exoplanet standards;
    % it is meant as a visibility/imaging stress test, not a realistic
    % detection scenario.

    p = struct();
    p.star_flux = 1.0;
    p.star_sigma_mas = 0.10;
    p.planet_sep_mas = 3.0;
    p.planet_pa_deg = 35.0;
    p.planet_sigma_mas = 0.08;
    p.contrast = 1e-1;

    targetCases.star_planet = p;


    %% Gaussian disk
    % Smooth circular resolved source.
    %
    % sigma_mas = 0.25 mas is deliberately comparable to the nominal
    % angular scale of a 1 Glambda array, so the visibility amplitude should
    % fall with rho_uv but not disappear immediately.

    p = struct();
    p.l0_mas = 0.0;
    p.m0_mas = 0.0;
    p.sigma_mas = 0.25;
    p.flux = 1.0;

    targetCases.gaussian_disk = p;


    %% Stellar disk
    % Circular stellar disk with linear limb darkening.
    %
    % The radius is chosen to be resolved but still compact relative to the
    % 40 mas half-FOV test image.

    p = struct();
    p.l0_mas = 0.0;
    p.m0_mas = 0.0;
    p.radius_mas = 0.45;
    p.limb_darkening = 0.4;
    p.flux = 1.0;

    targetCases.stellar_disk = p;


    %% Asteroid-like ellipse
    % Elongated resolved target with optional edge darkening and brightness
    % gradient.
    %
    % This is the most useful early NEO-style morphology. The dimensions are
    % large enough to be resolved by the current ~1 Glambda test coverage.

    p = struct();
    p.l0_mas = 0.0;
    p.m0_mas = 0.0;
    p.a_mas = 0.80;
    p.b_mas = 0.35;
    p.pa_deg = 25.0;
    p.flux = 1.0;
    p.limb_darkening = 0.35;
    p.gradient = 0.45;

    targetCases.ellipse = p;


    %% Crescent
    % Crescent-like extended target made by subtracting an offset inner disk
    % from an outer disk.
    %
    % This is intentionally larger than the Gaussian/ellipse tests because
    % crescent morphology needs enough pixels to show the shape.

    p = struct();
    p.l0_mas = 0.0;
    p.m0_mas = 0.0;
    p.outer_radius_mas = 1.20;
    p.inner_radius_mas = 0.95;
    p.inner_offset_l_mas = 0.42;
    p.inner_offset_m_mas = -0.08;
    p.flux = 1.0;
    p.smooth_sigma_pix = 1.2;

    targetCases.crescent = p;


    %% Multiple compact sources
    % Arbitrary collection of compact sources.
    %
    % These are separated by several mas so that the phase behavior is clear
    % and the dirty image should show multiple shifted beam responses.

    p = struct();

    p.sources(1).l_mas = -3.0;
    p.sources(1).m_mas =  1.0;
    p.sources(1).sigma_mas = 0.10;
    p.sources(1).flux = 1.0;

    p.sources(2).l_mas =  2.5;
    p.sources(2).m_mas = -1.5;
    p.sources(2).sigma_mas = 0.12;
    p.sources(2).flux = 0.6;

    p.sources(3).l_mas =  0.8;
    p.sources(3).m_mas =  3.5;
    p.sources(3).sigma_mas = 0.08;
    p.sources(3).flux = 0.25;

    targetCases.multi_point = p;
end