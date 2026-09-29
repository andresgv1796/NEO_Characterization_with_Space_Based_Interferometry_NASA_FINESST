function qpo = load_qpo_torus_from_family(qpoFamilyFile, qpoFamilyVariable, qpoIndex, N, sys)
%LOAD_QPO_TORUS_FROM_FAMILY Load one already-corrected GMOS family member.

    if ~isfile(qpoFamilyFile)
        located = which(qpoFamilyFile);
        if isempty(located)
            error('Could not find QPO family file "%s" on path.', qpoFamilyFile);
        end
        qpoFamilyFile = located;
    end

    S = load(qpoFamilyFile);
    if isfield(S, qpoFamilyVariable)
        fam = S.(qpoFamilyVariable);
    else
        names = fieldnames(S);
        if numel(names) ~= 1
            error('Could not infer family variable. Set qpoFamilyVariable explicitly.');
        end
        fam = S.(names{1});
        qpoFamilyVariable = names{1};
    end

    zQPO = get_family_member_z(fam, qpoIndex);

    if numel(zQPO) == 6*N + 3
        fprintf('Detected [z; mu] vector in family. Using mu from system object for analysis, stripping final entry.\n');
        zQPO = zQPO(1:end-1);
    end

    [XCurve, TQPO, rhoQPO] = gmos_unpack_z(zQPO, N);
    CCurve = Jacobi_Constant(XCurve.', sys.mu);

    qpo = struct();
    qpo.source = 'loadedFamily';
    qpo.file = qpoFamilyFile;
    qpo.variable = qpoFamilyVariable;
    qpo.index = qpoIndex;
    qpo.z = zQPO;
    qpo.XCurve = XCurve;
    qpo.T = TQPO;
    qpo.rho = rhoQPO;
    qpo.N = N;
    qpo.CCurve = CCurve;
    qpo.CMean = mean(CCurve);
    qpo.CMin = min(CCurve);
    qpo.CMax = max(CCurve);
    qpo.system = sys;
end
