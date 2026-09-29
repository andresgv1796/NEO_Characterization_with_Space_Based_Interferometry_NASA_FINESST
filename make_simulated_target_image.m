function simulated = make_simulated_target_image(simulatedType, nPix, fov_mas, p)
%MAKE_simulated_TARGET_IMAGE Create configurable simulated sky-brightness targets.
%
% simulated = make_simulated_target_image(simulatedType, nPix, fov_mas, p)
%
% Coordinates:
%   l, m are angular offsets from the target direction.
%   Units for all angular source parameters are milliarcseconds.
%   Position angles are measured counterclockwise from +l toward +m.
%
% Supported simulatedType:
%   'delta'
%   'point'
%   'binary'
%   'star_planet'
%   'gaussian_disk'
%   'ellipse'
%   'crescent'
%   'multi_point'
%   'stellar_disk'

    simulatedType = lower(char(simulatedType));

    mas2rad = pi/(180*3600*1000);

    l_mas = linspace(-fov_mas, fov_mas, nPix);
    m_mas = linspace(-fov_mas, fov_mas, nPix);

    [Lmas, Mmas] = meshgrid(l_mas, m_mas);

    I = zeros(nPix, nPix);

    switch simulatedType

        case 'delta'

            I(:) = 0;

            [~, i0] = min(abs(l_mas - p.l0_mas));
            [~, j0] = min(abs(m_mas - p.m0_mas));

            I(j0,i0) = p.flux;

        case 'point'
            I = add_gaussian_source( ...
                I, Lmas, Mmas, ...
                p.l0_mas, p.m0_mas, p.sigma_mas, p.flux);

        case 'binary'
            dx = 0.5*p.sep_mas*cosd(p.pa_deg);
            dy = 0.5*p.sep_mas*sind(p.pa_deg);

            l1 = p.center_l_mas - dx;
            m1 = p.center_m_mas - dy;

            l2 = p.center_l_mas + dx;
            m2 = p.center_m_mas + dy;

            I = add_gaussian_source( ...
                I, Lmas, Mmas, l1, m1, p.sigma1_mas, p.flux1);

            I = add_gaussian_source( ...
                I, Lmas, Mmas, l2, m2, p.sigma2_mas, p.flux2);

        case 'star_planet'
            lp = p.planet_sep_mas*cosd(p.planet_pa_deg);
            mp = p.planet_sep_mas*sind(p.planet_pa_deg);

            I = add_gaussian_source( ...
                I, Lmas, Mmas, ...
                0.0, 0.0, p.star_sigma_mas, p.star_flux);

            I = add_gaussian_source( ...
                I, Lmas, Mmas, ...
                lp, mp, p.planet_sigma_mas, p.star_flux*p.contrast);

        case 'gaussian_disk'
            X = Lmas - p.l0_mas;
            Y = Mmas - p.m0_mas;

            I = p.flux*exp(-(X.^2 + Y.^2)/(2*p.sigma_mas^2));

        case 'ellipse'
            X0 = Lmas - p.l0_mas;
            Y0 = Mmas - p.m0_mas;

            phi = deg2rad(p.pa_deg);

            X =  cos(phi)*X0 + sin(phi)*Y0;
            Y = -sin(phi)*X0 + cos(phi)*Y0;

            R = (X/p.a_mas).^2 + (Y/p.b_mas).^2;

            I(R <= 1) = p.flux;

            I = I .* max(0, 1 - p.limb_darkening*R);
            I = I .* (1 + p.gradient*X/p.a_mas);
            I(I < 0) = 0;

        case 'crescent'
            R_outer = sqrt((Lmas - p.l0_mas).^2 + ...
                           (Mmas - p.m0_mas).^2);

            R_inner = sqrt((Lmas - p.l0_mas - p.inner_offset_l_mas).^2 + ...
                           (Mmas - p.m0_mas - p.inner_offset_m_mas).^2);

            I(R_outer <= p.outer_radius_mas) = p.flux;
            I(R_inner <= p.inner_radius_mas) = 0;

            I = smooth_image_simple(I, p.smooth_sigma_pix);

        case 'multi_point'
            for k = 1:numel(p.sources)
                I = add_gaussian_source( ...
                    I, Lmas, Mmas, ...
                    p.sources(k).l_mas, ...
                    p.sources(k).m_mas, ...
                    p.sources(k).sigma_mas, ...
                    p.sources(k).flux);
            end

        case 'stellar_disk'
            
            X = Lmas - p.l0_mas;
            Y = Mmas - p.m0_mas;

            R = sqrt(X.^2 + Y.^2);
            rNorm = R/p.radius_mas;

            inside = rNorm <= 1;

            mu = zeros(size(R));
            mu(inside) = sqrt(1 - rNorm(inside).^2);

            I(inside) = p.flux*(1 - p.limb_darkening*(1 - mu(inside)));


       

        otherwise
            error('Unknown simulatedType "%s".', simulatedType);
    end

    I = I/sum(I(:));

    simulated = struct();
    simulated.type = simulatedType;
    simulated.params = p;
    simulated.I = I;

    simulated.l_mas = l_mas;
    simulated.m_mas = m_mas;
    simulated.L_mas = Lmas;
    simulated.M_mas = Mmas;

    simulated.l_rad = l_mas*mas2rad;
    simulated.m_rad = m_mas*mas2rad;
    simulated.L_rad = Lmas*mas2rad;
    simulated.M_rad = Mmas*mas2rad;

    simulated.fov_mas = fov_mas;
    simulated.nPix = nPix;
end

 