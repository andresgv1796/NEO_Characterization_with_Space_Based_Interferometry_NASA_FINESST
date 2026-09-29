function cleanMask = make_circular_clean_mask(l_mas, m_mas, l0_mas, m0_mas, radius_mas)
%MAKE_CIRCULAR_CLEAN_MASK
%
% Create a circular CLEAN mask centered at (l0_mas, m0_mas).

    [Lmas, Mmas] = meshgrid(l_mas, m_mas);

    R = sqrt((Lmas - l0_mas).^2 + (Mmas - m0_mas).^2);

    cleanMask = R <= radius_mas;
end