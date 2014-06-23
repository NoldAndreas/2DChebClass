function [V1,DV1] = getV1DV1(x1,t,optsPhys)

%--------------------------------------------------------------------------
% Set up coordinates
%--------------------------------------------------------------------------

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
else  % d=2
    x1=C(:,1);
    x2=C(:,2);
end
               
%--------------------------------------------------------------------------
% Get potential and gradient structures
%--------------------------------------------------------------------------

V1DV1=str2func(optsPhys.V1DV1);

if(nargout(V1DV1) == 3)
    confined = true;
else
    confined = false;
end

if(d==1)      % 1D potential
    if(confined)
        [VBackStruct,VAddStruct,VGeomStruct] = V1DV1(x1,t,optsPhys);
    else
        [VBackStruct,VAddStruct] = V1DV1(x1,t,optsPhys);
    end
else                  % 2D potential
    if(confined)
        [VBackStruct,VAddStruct,VGeomStruct] = V1DV1(x1,x2,t,optsPhys);
    else
        [VBackStruct,VAddStruct] = V1DV1(x1,x2,t,optsPhys);
    end
end

%--------------------------------------------------------------------------
% Extract potential
%--------------------------------------------------------------------------

VBack = VBackStruct.V;
VAdd  = VAddStruct.V;

% add to give total potential
V1= VBack + VAdd;

if(confined)
    VG = VGeomStruct.V;
    V1 = V1 + VG;
end

%--------------------------------------------------------------------------
% Construct gradient depending on calculation type, geometry and dimension
%--------------------------------------------------------------------------

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
            DV1 = DVBack + DVAdd;
            
            if(confined)
                DV1 = DV1 + (VGeomStruct.DV)';
            end
            
            % multiply by 1/R
            R = x1';
            DV1 = DV1./R;
            DV1(R==0) = 0;

            % replicate the rows dim times (dim,nParticles)
            DV1 = DV1(ones(1,dim),:);

        case 'planar'  % dim 1,2,3,  d = 1
            % gradient should act in potDir direction

            nParticles=optsPhys.nParticles;

            DV1 = zeros(dim,nParticles);
            DV1(end,:) = (VBackStruct.DV)' + (VAddStruct.DV)';
            
            if(confined)
                DV1(end,:) = DV1(end,:) + (VGeomStruct.DV)';
            end
            
        otherwise

            nParticles = optsPhys.nParticles;

            DV1 = zeros(dim,nParticles);

    end   % geometries

    % convert back to column vectors of size (nParticles x dim,1)
    DV1 = DV1(:);

    if(strcmp(geom,'spherical'))
        DV1 = DV1.*x;
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
            DV11 = VBackStruct.dy1 + VAddStruct.dy1;
            DV12 = VBackStruct.dy2 + VAddStruct.dy2;
            
            if(confined)
                DV11 = DV11 + VGeomStruct.dy1;
                DV12 = DV12 + VGeomStruct.dy2;
            end
            
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

            DV1R   = VBackStruct.dy1 + VAddStruct.dy1;
            DV1T_r = VBackStruct.dy2 + VAddStruct.dy2;
            
            if(confined)
                DV1R   = DV1R + VGeomStruct.dy1;
                DV1T_r = DV1T_r + VGeomStruct.dy2;
            end
            
            T = C(:,2);
            cosT = cos(T);
            sinT = sin(T);

            DV11 = DV1R .* cosT - DV1T_r .* sinT;
            DV12 = DV1R .* sinT + DV1T_r .* cosT;
            
    end    % geometry switch

    %--------------------------------------------------------------
    % Construct gradient from coordinate derivatives
    %--------------------------------------------------------------

    % put each derivative into a row (2,nParticles)
    DV1 = [DV11' ; DV12'];
    % convert to (nParticles,1) with correct x1, x2
    % ordering
    DV1 = DV1(:);

end    % dimension switch
      
end
