function dirty = make_dirty_image_from_observations(obs, l_rad, m_rad, addHermitian, weightMode)
%MAKE_DIRTY_IMAGE_FROM_OBSERVATIONS
%
% Form the dirty beam and dirty image from nonuniform visibility samples.
%
% Inputs:
%   obs          : output from make_simulated_interferometric_observations
%   l_rad        : image-grid l coordinates [rad]
%   m_rad        : image-grid m coordinates [rad]
%   addHermitian : true to append (-u,-v,V*) counterparts
%   weightMode   : visibility weighting mode
%                  'natural' : equal visibility weights
%                  'uniform' : inverse local UV-density weights
%                  'radial'  : weights proportional to rho_uv
%                  'radial2' : weights proportional to rho_uv^2
%
% Output:
%   dirty contains the dirty image, dirty beam, and effective UV set.
%
% Forward model:
%   V(u,v) = sum I(l,m) exp[-2*pi*i*(u*l + v*m)]
%
% Dirty-image synthesis:
%   I_D(l,m) = (1/sum w_k) sum w_k V_k exp[+2*pi*i*(u_k*l + v_k*m)]
%
% Dirty beam:
%   B_D(l,m) = (1/sum w_k) sum w_k exp[+2*pi*i*(u_k*l + v_k*m)]

    u = obs.u(:);
    v = obs.v(:);
    V = obs.V(:);

    if addHermitian
        u = [u; -u];
        v = [v; -v];
        V = [V; conj(V)];
    end

    w = compute_dirty_visibility_weights(u, v, l_rad, m_rad, weightMode);

    %% Enforce exact Hermitian weight symmetry

    if addHermitian

        nPhysical = ...
            numel(obs.u);

        if numel(w) ~= 2*nPhysical
            error([ ...
                'Hermitian weighting expected exactly twice ', ...
                'the physical visibility count.']);
        end

        wPhysical = ...
            w(1:nPhysical);

        wConjugate = ...
            w(nPhysical+1:end);

        wPair = ...
            0.5*(wPhysical + wConjugate);

        w = [ ...
            wPair; ...
            wPair];

    end
    sumWeights = sum(w);

    nL = numel(l_rad);
    nM = numel(m_rad);

    dirtyImage = zeros(nM, nL);
    dirtyBeam  = zeros(nM, nL);

    lCol = l_rad(:);

    parfor r = 1:nM

        mVal = m_rad(r);

        phase = 2*pi*(lCol*u.' + mVal*v.');
        E = exp(1i*phase);

        dirtyRow = E*(w .* V)/sumWeights;
        beamRow  = E*w/sumWeights;

        dirtyImage(r,:) = dirtyRow.';
        dirtyBeam(r,:)  = beamRow.';
    end

    [Lrad, Mrad] = meshgrid(l_rad, m_rad);

    dirty = struct();

    dirty.l_rad = l_rad;
    dirty.m_rad = m_rad;
    dirty.L_rad = Lrad;
    dirty.M_rad = Mrad;

    dirty.l_mas = l_rad*(180/pi)*3600*1000;
    dirty.m_mas = m_rad*(180/pi)*3600*1000;
    dirty.L_mas = Lrad*(180/pi)*3600*1000;
    dirty.M_mas = Mrad*(180/pi)*3600*1000;

   
    dirty.weightMode = weightMode;

    dirty.imageComplex = dirtyImage;
    dirty.beamComplex = dirtyBeam;
    
    if addHermitian

        beamImagMax = ...
            max(abs(imag(dirtyBeam(:))));

        imageImagMax = ...
            max(abs(imag(dirtyImage(:))));

        beamScale = ...
            max(1,max(abs(real(dirtyBeam(:)))));

        imageScale = ...
            max(1,max(abs(real(dirtyImage(:)))));

        if beamImagMax > 1e-10*beamScale
            error( ...
                ['Hermitian dirty beam is not numerically real. ', ...
                'max imaginary component = %.6e.'], ...
                beamImagMax);
        end

        if imageImagMax > 1e-10*imageScale
            error( ...
                ['Hermitian dirty image is not numerically real. ', ...
                'max imaginary component = %.6e.'], ...
                imageImagMax);
        end

    end
    dirty.image = real(dirtyImage);
    dirty.beam = real(dirtyBeam);

    dirty.imageAbs = abs(dirtyImage);
    dirty.beamAbs = abs(dirtyBeam);

    dirty.addHermitian = addHermitian;

    dirty.u = u(:);
    dirty.v = v(:);
    dirty.V = V(:);

    dirty.weights = w(:);
    dirty.sumWeights = sumWeights;
 

end


function w = compute_dirty_visibility_weights(u, v, l_rad, m_rad, weightMode)
%COMPUTE_DIRTY_VISIBILITY_WEIGHTS
%
% Visibility weighting for dirty imaging.
%
% Notes:
%   natural weighting gives every measured visibility equal weight.
%
%   uniform weighting approximates inverse local UV-density weighting by
%   binning samples on the Fourier grid associated with the image-plane
%   field of view. Dense UV cells receive lower weight.
%
%   radial and radial2 are diagnostic long-baseline-emphasizing weights.

    weightMode = char(weightMode);

    switch lower(weightMode)

        case 'natural'

            w = ones(size(u));

        case 'uniform'

            w = compute_uniform_uv_weights(u, v, l_rad, m_rad);

        case 'radial'

            rho = hypot(u, v);
            rhoScale = max(rho);
            w = rho/rhoScale;

        case 'radial2'

            rho = hypot(u, v);
            rhoScale = max(rho);
            w = (rho/rhoScale).^2;

        otherwise

            error('Unknown dirty-image weighting mode: %s', weightMode);
    end

    w(~isfinite(w)) = 0;

    if sum(w) <= 0
        error('Dirty-image weights sum to zero for weightMode = %s.', weightMode);
    end
end


function w = compute_uniform_uv_weights( ...
    u, v, l_rad, m_rad)

%COMPUTE_UNIFORM_UV_WEIGHTS
%
% Uniform UV weighting based on inverse occupancy of Fourier-grid cells.
%
% The UV-cell coordinates are defined relative to the Fourier origin,
% rather than relative to min(u),min(v). This preserves the symmetry
%
%       (u,v) <-> (-u,-v)
%
% required by Hermitian-completed visibility data.

    u = u(:);
    v = v(:);

    nL = numel(l_rad);
    nM = numel(m_rad);

    deltaL = ...
        abs(l_rad(2) - l_rad(1));

    deltaM = ...
        abs(m_rad(2) - m_rad(1));


    %% Natural Fourier-cell spacing associated with image grid

    du = ...
        1/(nL*deltaL);

    dv = ...
        1/(nM*deltaM);


    %% Origin-centered integer UV-cell coordinates
    %
    % This guarantees, apart from floating-point tie cases,
    %
    %       iu(-u) = -iu(u)
    %       iv(-v) = -iv(v).

    iu = ...
        round(u/du);

    iv = ...
        round(v/dv);


    %% Find UV-cell occupancy

    uvCell = [ ...
        iu, ...
        iv];

    [~,~,cellMembership] = ...
        unique( ...
            uvCell, ...
            'rows');

    occupancy = ...
        accumarray( ...
            cellMembership, ...
            1);


    %% Uniform weighting
    %
    % Every visibility receives inverse local sampling density.

    w = ...
        1 ./ occupancy(cellMembership);


    %% Normalize mean weight to unity

    w = ...
        w/mean(w);

end