function res = FGS(sim, covar)
% Fast Gaussian Simulation (FGS) generates Gaussian fields.
%   |res = FGS(sim, covar)|
%
%   INPUT:
%   |sim| struct with simulation settings
%   |sim.s| : size of grid as a vector, defines the dimension of the
%   grid.(default: [100 100])
%   |sim.n| : Number of simulation to perform. (default: 1)
%   |sim.tol| : Accuracy of the covariance map. For computational reason 
%   the covariance map is only accounted up to a range where the normalized
%   covariance value is tol.(default: 1e-6)
%   |sim.DisplayProgression| : Display a progression bar (default: false)
%   |sim.seed| : Control random number generation: 'default', 'shuffle' or 
%   a number. (default: 'shuffle')
% 
%   |covar| struct with covariance function/variogram settings. See
%   |covarIni.m| for documentation.
% 
%   OUTPUT:
%   |res| cell of size [sim.n x 1] with simulation field


% Validate input
validateattributes(sim,{'struct'},{})
if ~isfield(sim, 's'),sim.s=[100 100]; end
validateattributes(sim.s,{'numeric'},{'vector','integer','positive'})
if ~isfield(sim, 'n'),sim.n=1; end
validateattributes(sim.n,{'numeric'},{'integer','positive','scalar'})
if ~isfield(sim, 'tol'),sim.tol=1e-6; end
validateattributes(sim.tol,{'numeric'},{'positive','nonzero','finite'})
if ~isfield(sim, 'DisplayProgression'),sim.DisplayProgression=false; end
if(~usejava('jvm')),sim.DisplayProgression=false; end
% validateattributes(sim.DisplayProgression,{},{'scalar'})
if ~isfield(sim, 'seed'),sim.seed='shuffle'; end

% Control random number generation
rng(sim.seed)

% Initialization of the covariance structure. See |covarIni.m| for doc.
covar = covarIni(covar);

% Define grid X
x = cell(numel(sim.s),1);
for i=1:numel(sim.s)
    x{i} = 1:sim.s(i);
end  
[X{1:numel(sim.s)}] = ndgrid(x{:});
X = reshape(cat(numel(sim.s)+1,X{:}),[],numel(sim.s));

% Find the number of grids needed to reach a covariance equal to tolerance
val=fsolve(@(h) covar.g(h)-sim.tol,1,optimset('Display','off','TolFun',sim.tol/10));
% ring=floor(0.5+(val*max(covar.range))./sim.s);
d = max(abs(mrdivide(val*eye(numel(sim.s)),covar.cx)));
ring = floor(0.5+d./sim.s);


% Initiate covariance kernel/matrix
K = zeros(prod(sim.s),1);

% Display bar
if(sim.DisplayProgression)
    h = waitbar(0,'covariance map...');
end

% Periodicity of the covariance
for l=1:prod(2*ring+1)
    [position{1:numel(sim.s)}]=ind2sub(2*ring+1,l);
    position2=([position{:}] - ring-1).*sim.s+ceil(sim.s/2);
    K=K+covar.gxx0(X,position2); 
    if(sim.DisplayProgression)
        waitbar(l/prod(2*ring+1),h)
    end
end

% Reshape K in n-d
K = reshape(K,sim.s);

if(sim.DisplayProgression)
    waitbar(0,h,'simulation(s)')
end

% Generate res in a cell
res = cell(sim.n,1);
fftnK = fftn(K).^.5;
for k=1:sim.n
    res{k}=real(ifftn(fftn(covar.var*randn(sim.s)).*fftnK));
    if(sim.DisplayProgression)
        h = waitbar(k/sim.n);
    end
end

if(sim.DisplayProgression)
    close(h);
end

end
