function [x,p]=stochasticDynamics(f,x0,p0,optsPhys,optsStoc,plotPosMask)
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

% extract physical values for ease of use
dim=optsPhys.dim;
kBT=optsPhys.kBT;
D0=optsPhys.D0;
m=optsPhys.m;
M=diag(m);
HI=optsStoc.HI;
noise=optsStoc.noise;
g=optsPhys.gamma;

flow = optsPhys.flow;

% get HI interaction function
if(HI)
    makeHI=str2func(optsStoc.HIType);
end

% extract time step information
tMax=optsStoc.tMax;
tSteps=optsStoc.tSteps;

% calculate dt
dt=tMax/tSteps;

times=(0:tSteps)*dt;

% number of particles
nParticles=length(x0)/dim;

% type of stochastic calculation, either Langevin ('rv') or Ermak-McCammon
% ('r')
type=optsStoc.type;

switch type % Langevin or EM dynamics
    
    case 'rv'
    %----------------------------------------------------------------------
    % Full Langevin dynamics
    %----------------------------------------------------------------------
        
    % matrices to store the positions and momenta
    x=zeros(dim*nParticles,tSteps+1);
    p=x;

    % set initial conditions
    x(:,1)=x0;
    p(:,1)=p0;

    % Gamma and A for HI=0, or immediately overwritten for HI=1
    %Gamma=g*eye(dim*nParticles);
    %A=sqrt(kBT*m*g)*eye(dim*nParticles);
    Gamma = diag(g);
    A = sqrt(kBT)*diag(sqrt(m.*g));
    
    % perform evolution
    for t=2:tSteps+1

        % calculate matrices if HI are on
        if(HI)
            %Gamma=makeGamma(x(:,t-1),optsPhys);
            Gamma=makeHI(x(:,t-1),optsPhys);
            
            if(noise)
                %A=sqrt(kBT*m)*chol(Gamma,'lower');
                [A,pdDim]=chol(Gamma*M,'lower');
                
                if(pdDim>0)
                    fprintf(1,'Gamma not positive definite\n');
                    fprintf(1,[num2str(t) '\n']);
                    eig(Gamma*diag(m))
                    x(:,:)=NaN;
                    pause
                    break;
                end
                
                A=sqrt(kBT)*A;
            end
        end

        % gradient of potential depends on positions
        DV=getDV(x(:,t-1),times(t),optsPhys);
                
        % perform evolution
        x(:,t) = x(:,t-1) + dt*p(:,t-1)./m;
        p(:,t) = p(:,t-1) - dt*( Gamma*p(:,t-1) + DV ) + sqrt(2)*sqrt(dt)*A*f(:,t-1);
 
    end

    
    case 'r' 
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

        if(HI)
            %D=makeD(x(:,t-1),optsPhys);
            D=makeHI(x(:,t-1),optsPhys);
           
            if(noise)
                [A,pdDim]=chol(D,'lower');
                
                if(pdDim>0)
                    fprintf(1,'D not positive definite\n');
                    fprintf(1,[num2str(t) '\n']);
                    z=x(3:3:end,t-1);
                    min(z)
                    x(:,:)=NaN;
                    break;
                end
            end
        end
        
        % gradient of potential depends on positions
        DV=getDV(x(:,t-1),times(t),optsPhys);

        % perform evolution -- see Ermak McCammon '78
        % assuming d{D_ij}/d{r_j}=0, which it is for Rotne-Prager        
        
        x(:,t) = x(:,t-1) - dt*D*DV/kBT + sqrt(2) *sqrt(dt)*A*f(:,t-1);
        
        if(flow)
            U = getFlow(x(:,t-1),times(t),optsPhys);
            x(:,t) = x(:,t) + U*dt;
        end
        
    end
    
    
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