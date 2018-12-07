
function covar = covarIni(covar)
% covarIni intiates the covariance structure.
%   |covar = covarIni(covar)|
%
%   |covar| struct with covariance function/variogram settings
%   |covar.model| : Covariance model. Choose amongst: 'nugget', 'triangle',
%   'circular', 'spherical', 'cubic', 'exponential', 'gaussian', 'stable', 
%   'power', 'k-bessel', 'logarithmic', 'cauchy', 'hyperbolic', 
%   'cardinal sin', 'matheron' .(default: none) 
%   |covar.var| : variance of the field. (default: 1)
%   |sim.range| : Range of the covariance function. Vector of size 
%   equal to dimension
%   |sim.azimuth| : Azimuth of the covariance. Vector of size equal to
%   size(|sim.range|)-1

% Checking input variable
validateattributes(covar,{'struct'},{})
validateattributes(covar.model,{'char'},{})
if ~isfield(covar, 'var')||isempty(covar.var),covar.var=1; end
validateattributes(covar.var,{'numeric'},{'nonnegative','scalar'})
validateattributes(covar.range,{'numeric'},{'vector','positive'})
if ~isfield(covar,'azimuth')||isempty(covar.azimuth),covar.azimuth=0; end
validateattributes(covar.azimuth,{'numeric'},{'vector'})

% Defines the covariance function |covar.g|
switch covar.model
    case 'nugget'
        covar.g = @(h,r) h==0;
    case 'triangle'
        assert(numel(covar.range)==1,'only valid in 1D')
        covar.g = @(h) max(1-h,0);
    case 'cosine'
        assert(numel(covar.range)==1,'only valid in 1D')
        covar.g = @(h) cos(h);
    case 'circular'
        assert(numel(covar.range)<=2,'only valid in 1D and 2D')
        covar.g = @(h) 2/pi*(acos(min(h,1))-min(h,1).*sqrt(1-min(h,1).^2));
    case 'spherical'
        assert(numel(covar.range)<=3,'only valid in 1D, 2D and 3D')
        covar.g = @(h) 1-3/2*min(h,1)+1/2*min(h,1).^3;
    case 'cubic'
        covar.g = @(h) 1 - 7*min(h,1).^2 + 35/4*min(h,1).^3 - 7/2*min(h,1).^5 + 3/4*min(h,1).^7;
    case 'hyperbolic'
        covar.g = @(h) 1./(1+h);
        warning('Approx the integrale')
    case 'exponential'
        covar.g = @(h) exp(-h);
    case 'gaussian'
        covar.g = @(h) exp(-h.^2);
    case 'stable'
        assert(isfield(covar, 'alpha'),'alpha covar is not properly define')
        assert(covar.alpha>0 && covar.alpha<=2,'Alpha value not possible')
        covar.g = @(h) exp(-(h).^covar.alpha);
    case 'power'
        assert(isfield(covar, 'alpha'),'alpha covar is not properly define')
        assert(covar.alpha>0 && covar.alpha<2,'Alpha value not possible')
        covar.g = @(h) 1-h.^covar.alpha;
        warning('Approx the integrale')
    case 'k-bessel'
        assert(isfield(covar, 'alpha'),'alpha covar is not properly define')
        assert(covar.alpha>0,'Alpha needs to be positive')
        covar.g = @(h) 1/(2^(covar.alpha-1) * gamma(covar.alpha)) .* max(h,eps).^covar.alpha .* besselk(covar.alpha,max(h,eps));
    case 'logarithmic'
        covar.g = @(h) 1-log(h+1);
    case 'cauchy'
        assert(isfield(covar, 'alpha'),'alpha covar is not properly define')
        assert(covar.alpha>0,'Alpha value not possible')
        covar.g = @(h) (1+h.^2).^covar.alpha;
    case 'cardinal sine'
        assert(numel(covar.range)<=3,'only valid in 1D, 2D and 3D')
        covar.g = @(h) sin(max(eps,h))./max(eps,h);
    case 'matheron'
        covar.g = @(h) 1./(h);
    otherwise
        error('Variogram type not defined')
end

% Build the rotational matrix
rot=eye(numel(covar.range));
cang=cosd(covar.azimuth); sang=sind(covar.azimuth);
for i=1:min(numel(covar.azimuth),numel(covar.range)-1)
    rot = rot * ( padarray(padarray([cang(i)-1,-sang(i);sang(i),cang(i)-1],[i-1 i-1],'pre'),[numel(covar.range)-1-i numel(covar.range)-1-i],'post')+eye(numel(covar.range)));
end

% Combine rotation and range (normalization)
covar.cx = rot/diag(covar.range);

% Defines the un-scaled covariance function calculation for
% covariance vector x vector (cross-covariance)
covar.gx = @(x) covar.g(squareform(pdist(x*covar.cx)))*covar.var;
% covariance vector x pt
covar.gxx0 = @(x,x0) covar.g(sqrt(sum((bsxfun(@minus,x,x0)*covar.cx).^2,2)))*covar.var;
% covariance vector1 x vector2
covar.gx1x2 = @(x1,x2) covar.g(pdist2(x1*covar.cx,x2*covar.cx))*covar.var;
end