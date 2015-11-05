function U = getFlow(x1,t,optsPhys)

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

UDU = str2func(optsPhys.UDU);


if(d==1)      % 1D potential
    disp('not implemented in 1D yet');
    pause    
else                  % 2D potential
    
    UStruct = UDU(x1,x2,t,optsPhys);
        
    U1 = UStruct.U1;
    U2 = UStruct.U2;

    % put each derivative into a row (2,nParticles)
    U = [U1' ; U2'];
    % convert to (nParticles,1) with correct x1, x2
    % ordering
    U = U(:);
end

end
