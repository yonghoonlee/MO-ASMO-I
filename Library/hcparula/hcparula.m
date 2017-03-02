function rgb = hcparula(varargin)
% High Contrast Parula-like Colormap Generator
% Author: Yong Hoon Lee (http://yonghoonlee.com, <a
% href="mailto:ylee196@illinois.edu">ylee196@illinois.edu</a>)
%
% This is a monotonic and linear grayscale intensity colormap with an
% increased contrast for grayscale printout readability.
% The colormap looks similar to Parula, but this code does not contain any
% formula, code, colormap scheme, or specific knowledge that MathWorks,
% Inc. generated or developed. Note that the original Parula colormap was
% created by MathWorks, Inc. and a default colormap in MATLAB.
%
%   RGB = MODIFIED_PARULA() generates a N=64 level colormap in RGB format.
%   RGB = MODIFIED_PARULA(N) generates a given N-level colormap in RGB
%   format.
%
%   See also COLORMAP, COLORBAR, RGBPLOT

    n = 64;
    low = 0.025;
    high = 0.925;
    if (nargin > 0)
        n = varargin{1};
    end
    xq = linspace(0,1,n);
    iq = linspace(low,high,n);
    rx = [0.00, 0.05, 0.10, 0.42, 0.82, 1.00];
    rf = [0.00, 0.08, 0.00, 0.00, 1.00, 1.00];
    rq = interp1(rx,rf,xq);
    gx = [0.00, 0.05, 0.42, 0.82, 0.96, 1.00];
    gf = [0.00, 0.00, 0.61, 0.77, 1.00, 1.00];
    gq = interp1(gx,gf,xq);
    bq = (iq - 0.2989*rq - 0.5870*gq)/0.1140;
    rgb = [rq',gq',bq'];
end