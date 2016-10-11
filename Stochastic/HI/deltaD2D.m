function D=deltaD2D(r,sigmaH)

% See von Hansen, Hinczewski, Netz, JCP, 134, 235102 (2011)

%--------------------------------------------------------------------------
% Extract x and z components 
%--------------------------------------------------------------------------

x = r(1:2:end);
z = r(2:2:end);

%--------------------------------------------------------------------------
% Calculate frequently-used terms
%--------------------------------------------------------------------------

nParticles = length(r)/2;

zp   = z';
z1z2 = z*zp;

nvec  = 1:nParticles;                    % row vector of column indices
nmask = nvec(ones(nParticles,1),:);     % nParticle rows of 1:nParticles (row length)

Rx = x(nmask)'-x(nmask);
Rz = z(nmask)'+z(nmask);  % note + as we're using reflected z

% equivalent to
% z1z2=bsxfun(@times,z,z');
% x1mx2=bsxfun(@minus,x,x');
% y1my2=bsxfun(@minus,y,y');
% z1pz2=bsxfun(@plus,z,z');

% 1/|r_i-r_j|
R1 = 1./( Rx.^2  + Rz.^2).^(1/2);
R3 = R1.*R1.*R1;   % multiplication faster than powers
R5 = R3.*R1.*R1;
R7 = R5.*R1.*R1;

z2 = zp(ones(nParticles,1),:); % z2 should be constant in columns
% equivalent to z2=repmat(z',nParticles,1); but faster

%--------------------------------------------------------------------------
% Calculate HI terms due to wall
%--------------------------------------------------------------------------

% Oseen of reflection (factor -1/(8 pi eta))
DORxx = R1 + Rx.^2.*R3;
DORzz = R1 + Rz.^2.*R3;
DORxz = Rx.*Rz.*R3;
DORzx = DORxz;

% Oseen wall terms (factor 1/(4 pi eta))
deltaDOxx = - z1z2.*( R3 - 3*Rx.^2.*R5 );
deltaDOzz =   z1z2.*( R3 - 3*Rz.^2.*R5 );
deltaDOxz = Rx.*( z2.*R3 - 3*z1z2.*Rz.*R5 );
deltaDOzx = Rx.*( z2.*R3 + 3*z1z2.*Rz.*R5 );

% Rotne-Prager of reflection (factor -a^2/(4 pi eta))
DRPRxx = 1/3*R3 - Rx.^2.*R5;
DRPRzz = 1/3*R3 - Rz.^2.*R5;
DRPRxz = - Rx.Rz.*R5;
DRPRzx = DRPRxz;

% Rotne-Prager wall terms (factor a^2/(4 pi eta))
deltaDRPxx = Rz.^2.*( R5   - 5*Rx.^2.*R7 );
deltaDRPzz = Rz.^2.*( 3*R5 - 5*Rz.^2.*R7 );
deltaDRPxz = - Rx.*Rz.*( 2*R5 - 5*Rz.^2.*R7 );
deltaDRPzx = - 5*Rx.*Rz.^3.*R7;

%--------------------------------------------------------------------------
% Fill appropriate matrix entries via masks
%--------------------------------------------------------------------------

c4 = 3/4*diag(sigmaH);
c8 = 3/8*diag(sigmaH);
ca24 = 3/16*diag(sigmaH).^3;

D=zeros(2*nParticles,2*nParticles);

% mask to reproduce 2x2 matrix nParticles x nParticles times
dvec  = (1:2)';                          % column vector of 1:dim
dmask = dvec(:,ones(nParticles,1));     % nParticle columns each 1:dim

xxMask = logical([1 0 ; 0 0]);
xxMask = xxMask(dmask,dmask);

zzMask = logical([0 0; 0 1]);
zzMask = zzMask(dmask,dmask);

xzMask = logical([0 1; 0 0]);
xzMask = xzMask(dmask,dmask);

zxMask=logical([0 0; 1 0]);
zxMask=zxMask(dmask,dmask);

D(xxMask) = -c8.*DORxx - ca24.*DRPRxx + c4.*deltaDOxx + ca24.*deltaDRPxx;

D(zzMask) = -c8.*DORzz - ca24.*DRPRzz + c4.*deltaDOzz + ca24.*deltaDRPzz;

D(xzMask) = -c8.*DORxz - ca24.*DRPRxz + c4.*deltaDOxz + ca24.*deltaDRPxz;

D(zxMask) = -c8.*DORzx - ca24.*DRPRzx + c4.*deltaDOzx + ca24.*deltaDRPzx;
