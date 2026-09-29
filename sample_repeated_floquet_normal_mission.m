function mission = sample_repeated_floquet_normal_mission( ...
    nSynodicOnePeriod,TPrimary,durationTU,nSamples,synodicPhase0_rad)
%SAMPLE_REPEATED_FLOQUET_NORMAL_MISSION Sample a repeated normal template.
%
% The one-period synodic Floquet-normal line is treated as the maintained
% reduced-order geometry. Deterministic low-discrepancy times sample a common
% physical mission horizon without re-integrating the STM or allocating every
% repeated primary-period node. The continuous synodic-to-inertial rotation is
% retained at every mission sample.

    narginchk(5,5);
    validateattributes(nSynodicOnePeriod,{'numeric'}, ...
        {'real','finite','2d','nrows',3},mfilename,'nSynodicOnePeriod',1);
    validateattributes(TPrimary,{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,'TPrimary',2);
    validateattributes(durationTU,{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,'durationTU',3);
    validateattributes(nSamples,{'numeric'}, ...
        {'real','finite','scalar','integer','>=',100},mfilename,'nSamples',4);
    validateattributes(synodicPhase0_rad,{'numeric'}, ...
        {'real','finite','scalar'},mfilename,'synodicPhase0_rad',5);

    nTemplate = size(nSynodicOnePeriod,2);
    if nTemplate < 20
        error('FINESST:MissionAverage:InsufficientTemplateSamples', ...
            'The one-period normal template requires at least 20 samples.');
    end
    normalNormError = max(abs(vecnorm(nSynodicOnePeriod,2,1)-1));
    if normalNormError > 1e-10
        error('FINESST:MissionAverage:NonunitTemplate', ...
            'The template-normal norm error reached %.3e.',normalNormError);
    end

    % A one-dimensional Kronecker sequence avoids resonance between a uniform
    % mission grid, TPrimary, and the rotating-frame period. The sequence is
    % deterministic and uniformly samples the requested physical time span.
    goldenConjugate = (sqrt(5)-1)/2;
    u = mod(((0:nSamples-1)+0.5)*goldenConjugate,1);
    t = durationTU*u;

    templatePhase = mod(t,TPrimary)/TPrimary;
    templateIndex = floor(nTemplate*templatePhase)+1;
    templateIndex = min(nTemplate,max(1,templateIndex));
    nSynodic = nSynodicOnePeriod(:,templateIndex);

    theta = synodicPhase0_rad+t;
    c = cos(theta);
    s = sin(theta);
    nInertial = [ ...
        c.*nSynodic(1,:)-s.*nSynodic(2,:); ...
        s.*nSynodic(1,:)+c.*nSynodic(2,:); ...
        nSynodic(3,:)];

    normError = max(abs(vecnorm(nInertial,2,1)-1));
    if normError > 1e-10
        error('FINESST:MissionAverage:UnitNormalFailure', ...
            'The sampled inertial-normal norm error reached %.3e.',normError);
    end

    mission = struct();
    mission.t = t;
    mission.nInertial = nInertial;
    mission.durationTU = durationTU;
    mission.nSamples = nSamples;
    mission.method = [ ...
        'deterministic Kronecker time quadrature of repeated one-period ' ...
        'synodic Floquet-normal template'];
end
