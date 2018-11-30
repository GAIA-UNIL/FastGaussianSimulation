
function covar = covarinitiate(covar)

validateattributes(covar,{'struct'},{})
validateattributes(covar.model,{'char'},{})
validateattributes(covar.c0,{'numeric'},{'positive','scalar'})
validateattributes(covar.range0,{'numeric'},{'vector','positive'})
validateattributes(covar.azimuth,{'numeric'},{'vector'})
assert(numel(covar.range0)==numel(covar.azimuth)+1,'Azimuth size need to be range-1')

switch covar.model
    case 'nugget'
        covar.g = @(h,r) h==0;
        intvario=1;
    case 'triangle'
        assert(numel(covar.range)==1,'only valid in 1D')
        intvario=1;
        covar.g = @(h) max(1-h,0);
    case 'circular'
        intvario=1.17;
        covar.g = @(h) 2/pi*(acos(min(h,1))-min(h,1).*sqrt(1-min(h,1).^2));
    case 'spherical'
        intvario=1.3;
        covar.g = @(h) 1-3/2*min(h,1)+1/2*min(h,1).^3;
    case 'cubic'
        intvario=1.43;
        covar.g = @(h) 1 - 7*min(h,1).^2 + 35/4*min(h,1).^3 - 7/2*min(h,1).^5 + 3/4*min(h,1).^7;
    case 'exponential'
        intvario = .41;
        covar.g = @(h) exp(-h);
    case 'gaussian'
        intvario=.58;
        covar.g = @(h) exp(-h.^2);
    case 'stable'
        assert(isfield(covar, 'alpha'),'alpha covar is not properly define')
        assert(covar.alpha>0 && covar.alpha<=2,'Alpha value not possible')
        intvario=.41;
        covar.g = @(h) exp(-(h).^covar.alpha);
    case 'power'
        assert(isfield(covar, 'alpha'),'alpha covar is not properly define')
        assert(covar.alpha>0 && covar.alpha<2,'Alpha value not possible')
        covar.g = @(h) 1-h.^covar.alpha;
        warning('Approx the integrale')
        intvario=1.5;
    case 'k-bessel'
        assert(isfield(covar, 'alpha'),'alpha covar is not properly define')
        %assert(covar.alpha>0 && covar.alpha<2,'Alpha value not possible')
        intvario=[.35 .5 .5];
        covar.g = @(h) 1/(2^(covar.alpha-1) * gamma(covar.alpha)) .* max(h,eps).^covar.alpha .* besselk(covar.alpha,max(h,eps));
    case 'logarithmic'
        covar.g = @(h) 1-log(h+1);
        warning('Approx the integrale')
        intvario=.7;
    case 'cauchy'
        assert(isfield(covar, 'alpha'),'alpha covar is not properly define')
        assert(covar.alpha>0,'Alpha value not possible')
        covar.g = @(h) (1+h.^2).^covar.alpha;
        warning('Approx the integrale')
        intvario=1;
    case 'hyperbolic'
        covar.g = @(h) 1./(1+h);
        warning('Approx the integrale')
        intvario=[.2 .05];
    case 'cardinal sine'
        intvario=.2;
        covar.g = @(h) sin(max(eps,h))./max(eps,h);
    case 'matheron'
        covar.g = @(h) 1./(h);
    otherwise
        error('Variogram type not defined')
end

if numel(covar.range0)==1 || numel(intvario)==1
    covar.range=covar.range0*intvario(1);
elseif numel(covar.range0)==2
    covar.range=covar.range0*intvario(2);
elseif numel(covar.range0)==3
    covar.range=covar.range0*intvario(3);
end

if numel(covar.range)==1 || numel(covar(1).azimuth)==0
    covar.cx = 1/diag(covar.range(1));
elseif numel(covar.range)==2
    cang=cosd(covar.azimuth); sang=sind(covar.azimuth);
    rot = [cang,-sang;sang,cang];
    covar.cx = rot/diag(fliplr(covar.range));
elseif numel(covar.range)==3
    cang=cosd(covar.azimuth); sang=sind(covar.azimuth);
    rot = [cang(1) -sang(1) 0;sang(1) cang(1) 0; 0 0 1] * [cang(2) 0 -sang(2);0 1 0;sang(2) 0 cang(2)];
    covar.cx = rot/diag(fliplr(covar.range));
end

covar.gx = @(x) covar.g(squareform(pdist(x*covar.cx)))*covar.c0;
covar.gxx0 = @(x,x0) covar.g(sqrt(sum((bsxfun(@minus,x,x0)*covar.cx).^2,2)))*covar.c0;
covar.gx1x2 = @(x1,x2) covar.g(pdist2(x1*covar.cx,x2*covar.cx))*covar.c0;
end