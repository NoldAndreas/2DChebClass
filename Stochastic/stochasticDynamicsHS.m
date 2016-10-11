function [x,p]=stochasticDynamicsHS(f,x0,p0,optsPhys,optsStoc,plotPosMask)
% [x,p]=Dynamics(f,x0,p0,optsPhys,optsStoc,plotPosMask)
%   performs stochastic dynamics calculattion with noise f and initial
%   positions x0 and momenta p0
%
% INPUTS: 
%  f        -- matrix of size dim*nParticles x tSteps giving the noise
%  x        -- initial positions, vector length dim*nParticles
%  p        -- initial momenta, vector length dim*nParticles
%  optsPhys -- a structure of size [1 1] containing at least:
%         dim           (dimension; 1 or 3)
%         D0            (diffusion constant)
%         kBT           (temperature)
%         nParticles    (number of particles)
%         sigma         (particle diameter)
%         m             (particle mass)
%         HI            (hydrodynamic interactions true/false)
%         V1DV1         (string identifying potential and gradient
%                        function)
%         [pot params]  (parameters for potential)
%
% optsStoc -- a structure of size [1 1] containing at least
%         tMax          (max time)
%         tSteps        (number of time steps)
%         type          ('rv'/'r')
%
% OUTPUTS:
%  x           -- a matrix of positions of size 
%                 (dim x nParticles, length(plotTimes))
%  p           -- a matrix of momenta of size 
%                 (dim x nParticles, length(plotTimes))

% remove HS part of potential - dealt with via resampling collisions
optsPhys.V2DV2 = 'free';

% extract physical values for ease of use
dim=optsPhys.dim;
kBT=optsPhys.kBT;
D0=optsPhys.D0;
sigma = optsPhys.sigma;
g=optsPhys.gamma;

% extract time step information
tMax=optsStoc.tMax;
tSteps=optsStoc.tSteps;

% calculate dt
dt=tMax/tSteps;

times=(0:tSteps)*dt;

% number of particles
nParticles=length(x0)/dim;

%----------------------------------------------------------------------
% Ermak-McCammon dynamics
%----------------------------------------------------------------------

% matrices to store the positions and momenta
x=zeros(dim*nParticles,tSteps+1);
p=x;

% set initial conditions
x(:,1)=x0;

% D and A for HI=0
%D=D0*eye(dim*nParticles);
%A=sqrt(D0)*eye(dim*nParticles);
D=diag(D0);
A=diag(sqrt(D0));

% perform evolution
for t=2:tSteps+1

    % gradient of potential depends on positions
    DV=getDV(x(:,t-1),times(t),optsPhys);

    v = - dt*D*DV/kBT;
    
    dr = v +  sqrt(2*dt)*A*f(:,t-1);
    
    drs = dr;
    
    % perform evolution -- see Ermak McCammon '78
    r = x(:,t-1) + dr;

    pairs = findOverlaps(r,sigma,dim);
    
    nOverlaps = size(pairs,1);
    
    for iOverlap = 1:nOverlaps
        p1 = pairs(iOverlap,1);
        p2 = pairs(iOverlap,2);
        
        pos1 = dim*(p1-1)+1:dim*(p1-1)+dim;
        pos2 = dim*(p2-1)+1:dim*(p2-1)+dim;
        
        r1 = x(pos1,t-1);
        r2 = x(pos2,t-1);
        
        [r12,n12] = getRij(r1,r2,dim);
        n12 = squeeze(n12);
        e12 = - n12;
        
        q0 = r12 - sigma(p1,p2);
        
        D1 = D(dim*p1,dim*p1);
        D2 = D(dim*p2,dim*p2);
        
        Dq = D1 + D2;
        
        v12 = v(pos1) - v(pos2);
        v12 = sum(v12.*n12);
        [q,~] = HSCollisionSampling(dt,q0,v12,Dq);
        
        dq = q - q0;
        
        dr1 = dr(pos1);
        dr2 = dr(pos2);
        
        dr21e = sum( (dr2 - dr1) .* e12 );
        
        dr1s = dr1 + D1/Dq * ( dr21e - dq) * e12;
        dr2s = dr2 - D2/Dq * ( dr21e - dq) * e12;
        
        drs(pos1) = dr1s;
        drs(pos2) = dr2s;
            
    end
    
    x(:,t) = x(:,t-1) + drs;
    
%     if(nOverlaps>0)
%         
%         pairsNew = findOverlaps(x(:,t),sigma,dim);
%         
%         pairs;
%         
%         if(~isempty(pairsNew))
%             disp('overlap!');
%         else
%             disp('ok');
%         end
%           
%     end

    
end
    
% if we don't specify any plotting/save times, do it at all of them.  Not
% really recommended as the files get huge and you can also easily run out
% of memory
if (isempty(plotPosMask))
    plotPosMask=true(tSteps+1,1);
end

% return only the values at the requested (plotPosMask) times
x=x(:,plotPosMask);
p=p(:,plotPosMask);