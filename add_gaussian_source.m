function I = add_gaussian_source(I, Lmas, Mmas, l0_mas, m0_mas, sigma_mas, flux)

    G = exp(-((Lmas - l0_mas).^2 + (Mmas - m0_mas).^2)/(2*sigma_mas^2));
    G = G/sum(G(:));

    I = I + flux*G;
end

