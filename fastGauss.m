function res = fastGauss(sim, covar)
% FastGauss Generate Gaussian fields.
%   |res = FastGauss(sim, covar)|
%
%   |sim| struct with simulation settings
%   |sim.s| : size of grid as a vector, defines the dimension of the
%   grid.(default: [100 100])
%   |sim.n| : Number of simulation to perform. (default: 1)
%   |sim.tol| : tolerance of the covariance map to include (see doc for
%   details). (default: 0.1)
% 
%   |covar| struct with covariance function/variogram settings. See
%   |covarinitiate| for documentation.
% 
%   |res| cell of size [sim.n x 1] with simulation field


% Validate input
validateattributes(sim,{'struct'},{})
if ~isfield(sim, 's'),sim.s=[100 100]; end
validateattributes(sim.s,{'numeric'},{'vector','integer','positive'})
if ~isfield(sim, 'n'),sim.n=1; end
validateattributes(sim.n,{'numeric'},{'integer','positive','scalar'})
if ~isfield(sim, 'tol'),sim.tol=0.1; end
validateattributes(sim.tol,{'numeric'},{'positive','nonzero','finite'})

covar = covarinitiate(covar);

% Define grid 
x = cells(numel(sim.s),1);
for i=1:numel(sim.s)
    x{i} = 1:sim.s(i);
end  
[X{1:numel(sim.s)}] = ndgrid(x{:});
X = reshape(permute([X{:}],[1 3 2]),[],numel(sim.s));

% Find the number of grids needed to reach a covariance equal to tolerance
val=fsolve(@(h) covar.g(h)-sim.tol,1,optimset('Display','off'));
ring=floor((val*max(covar.range))./sim.s);

% Initiate covariance kernel/matrix
K = zeros(prod(sim.s),1);

% Periodicity of the covariance
for l=1:prod(2*ring+1)
    [position{1:numel(sim.s)}]=ind2sub(2*ring+1,l);
    position2=([position{:}] - ring-1).*sim.s+ceil(sim.s/2);
    K=K+covar.gxx0(X,position2); 
end

% Reshape K in n-d
K = reshape(K,sim.s);

% Generate res in a cell
res = cell(sim.n,1);
fftnK = fftn(K).^.5;
for k=1:sim.n
    res{k}=real(ifftn(fftn(randn(sim.s)).*fftnK));
end

end
