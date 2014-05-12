function [V1,DV1,VBack,DVBack,VAdd]=getV1DV1(x1,x2,t,optsPhys)

%--------------------------------------------------------------------------
% Set up coordinates
%--------------------------------------------------------------------------

switch optsPhys.type
    
    case 'stoc'
        
        % ensure column vector, e.g. for slice sampling which uses row
        % vectors
        x=x1(:);

        dim=optsPhys.dim;     % stochastic dimension (1, 2 or 3)
        geom=optsPhys.geom;   % geometry of potential
        
        % get coordinates for each particle (of size nParticles x d)
        C = getCoords(x,dim,geom);

        d = size(C,2);        % geometry dimension
        
        if (d==1)
            x = x1;
            x1=C;
            x2=[];
        else  % d=2
            x1=C(:,1);
            x2=C(:,2);
        end
               
    otherwise
        
        % input is already in coordinate form
        
end

%--------------------------------------------------------------------------
% Get potential and gradient structures
%--------------------------------------------------------------------------

V1DV1=str2func(optsPhys.V1DV1);

if(isempty(x2))       % 1D potential
    [VBackStruct,VAddStruct]=V1DV1(x1,t,optsPhys);
else                  % 2D potential
    [VBackStruct,VAddStruct]=V1DV1(x1,x2,t,optsPhys);
end
    
    
%--------------------------------------------------------------------------
% Extract potential
%--------------------------------------------------------------------------

VBack = VBackStruct.V;
VAdd  = VAddStruct.V;

% add to give total potential
V1= VBack + VAdd;

%--------------------------------------------------------------------------
% Construct gradient depending on calculation type, geometry and dimension
%--------------------------------------------------------------------------

switch optsPhys.type
    
    %----------------------------------------------------------------------
    % Stochastic calculations
    %----------------------------------------------------------------------
    
    case 'stoc'
        
        % potential should just be a number -- only used in slice sampling
        % the initial/final distributions
        
        % sum over all particles
        V1=sum(V1);

        %------------------------------------------------------------------
        % 1-dimensional potentials
        %------------------------------------------------------------------
        
        if (d==1)
            
            switch geom
                
                %----------------------------------------------------------
                % Spherical
                %----------------------------------------------------------
                
                case 'spherical'  % dim 3,  d = 1
                  
                    % get gradients as row vectors (1,nParticles)
                    DVBack = (VBackStruct.DV)';  
                    DVAdd  = (VAddStruct.DV)';
                    
                    % multiply by 1/R
                    R = x1';
                    DVBack = DVBack./R;
                    DVAdd  = DVAdd./R;
                    DVBack(R==0) = 0;
                    DVAdd(R==0) = 0;
                    
                    % replicate the rows dim times (dim,nParticles)
                    DVBack = DVBack(ones(1,dim),:);
                    DVAdd  = DVAdd(ones(1,dim),:);

                case 'planar'  % dim 1,2,3,  d = 1
                    % gradient should act in potDir direction
                    
                    nParticles=optsPhys.nParticles;
                    
                    DVBack = zeros(dim,nParticles);
                    DVAdd  = zeros(dim,nParticles);

                    DVBack(end,:) = (VBackStruct.DV)';
                    DVAdd(end,:)  = (VAddStruct.DV)';

%                   ADD IN CASES FOR DIFFERENT PLANAR DIRECTIONS?           
%                     if(dim==1)
%                         potDir=1;
%                     else
%                         potDir=optsPhys.potDir;
%                     end
%                     
%                     DVBack(potDir,:) = (VBackStruct.grad)';
%                     DVAdd(potDir,:)  = (VAddStruct.grad)';
                    
                otherwise
                    
                    nParticles=optsPhys.nParticles;
                    
                    DVBack=zeros(dim,nParticles);
                    DVAdd=zeros(dim,nParticles);
                        
            end   % geometries
                    
            % convert back to column vectors of size (nParticles x dim,1)
            DVBack = DVBack(:);
            DVAdd  = DVAdd(:);
            
            if(strcmp(geom,'spherical'))
                DVBack = DVBack.*x;
                DVAdd  = DVAdd.*x;
            end


        %------------------------------------------------------------------
        % 2-dimensional potentials
        %------------------------------------------------------------------
            
        else % d=2
            
            switch geom
                
                %----------------------------------------------------------
                % Planar 2D
                %----------------------------------------------------------
                
                case 'planar2D'  % dim = d = 2
                    DVBack1 = VBackStruct.dy1;
                    DVBack2 = VBackStruct.dy2;

                    DVAdd1  = VAddStruct.dy1;
                    DVAdd2  = VAddStruct.dy2;
                    
                %----------------------------------------------------------
                % Polar 2D
                %----------------------------------------------------------    
                    
                case 'polar2D'  % dim = d = 2
                    % need to convert from polar gradient to cartesian
                    % gradient
                    
                    % note that the theta derivative obtained from V1DV1
                    % should actually be 1/r d_theta
                     
                    % x = r cos t, y = r sin t
                    % dV/dx1 = dV/dr cos t - (dV/dt sin t) / r
                    % dV/dx2 = dV/dr sin t + (dV/dt cos t) / r
                    
                    % Coordinates are in the order r, theta
                    
                    DVBackR   = VBackStruct.dy1;
                    DVBackT_r = VBackStruct.dy2;

                    DVAddR    = VAddStruct.dy1;
                    DVAddT_r  = VAddStruct.dy2;
                    
                    T = C(:,2);
                    cosT = cos(T);
                    sinT = sin(T);
                    
                    DVBack1 = DVBackR .* cosT - DVBackT_r .* sinT;
                    DVBack2 = DVBackR .* sinT + DVBackT_r .* cosT;
                    
                    DVAdd1  = DVAddR .* cosT  - DVAddT_r .* sinT;
                    DVAdd2  = DVAddR .* sinT  + DVAddT_r .* cosT;
                    
            end    % geometry switch
            
            %--------------------------------------------------------------
            % Construct gradient from coordinate derivatives
            %--------------------------------------------------------------
            
            % put each derivative into a row (2,nParticles)
            DVBack = [DVBack1' ; DVBack2'];
            % convert to (nParticles,1) with correct x1, x2
            % ordering
            DVBack = DVBack(:);

            % do same for DVAdd
            DVAdd  = [DVAdd1' ; DVAdd2'];
            DVAdd  = DVAdd(:);
                    
        end    % dimension switch
    
    %----------------------------------------------------------------------
    % DDFT calculations
    %----------------------------------------------------------------------

    otherwise  % eg DDFT and plotting
        
        DVBack = VBackStruct.DV;
        DVAdd  = VAddStruct.DV;
    
end   % calculation type


% add to give total gradient
DV1=DVBack + DVAdd;
   
end
