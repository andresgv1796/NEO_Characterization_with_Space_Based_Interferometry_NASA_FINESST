function Iout = smooth_image_simple(I, sigmaPix)
%SMOOTH_IMAGE_SIMPLE Apply a separable Gaussian blur.
%
% Used only to soften pixel-sharp toy target edges.

    halfWidth = ceil(3*sigmaPix);
    x = -halfWidth:halfWidth;

    g = exp(-0.5*(x/sigmaPix).^2);
    g = g/sum(g);

    Iout = conv2(conv2(I, g, 'same'), g.', 'same');
end