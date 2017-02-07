function D=deltaD3D(r,sigmaH)

% See von Hansen, Hinczewski, Netz, JCP, 134, 235102 (2011)

%--------------------------------------------------------------------------
% Extract x and z components 
%--------------------------------------------------------------------------

x = r(1:3:end);
y = r(2:3:end);
z = r(3:3:end);

%--------------------------------------------------------------------------
% Calculate frequently-used terms
%--------------------------------------------------------------------------

nParticles = length(r)/4;

zp   = z';
z1z2 = z*zp;

nvec  = 1:nParticles;                    % row vector of column indices
nmask = nvec(ones(nParticles,1),:);     % nParticle rows of 1:nParticles (row length)

Rx = x(nmask)'-x(nmask);
Ry = y(nmask)'-y(nmask);
Rz = z(nmask)'+z(nmask);  % note + as we're using reflected z

% equivalent to
% z1z2=bsxfun(@times,z,z');
% x1mx2=bsxfun(@minus,x,x');
% y1my2=bsxfun(@minus,y,y');
% z1pz2=bsxfun(@plus,z,z');

% 1/|r_i-r_j|
R1 = 1./( Rx.^2  + Ry.^2 + Rz.^2).^(1/2);
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
DORyy = R1 + Ry.^2.*R3;
DORzz = R1 + Rz.^2.*R3;
DORxy = Rx.*Ry.*R3;
DORyx = DORxy;
DORxz = Rx.*Rz.*R3;
DORzx = DORxz;
DORyz = Ry.*Rz.*R3;
DORzy = DORyz;

% Oseen wall terms (factor 1/(4 pi eta))
deltaDOxx = - z1z2.*( R3 - 3*Rx.^2.*R5 );
deltaDOyy = - z1z2.*( R3 - 3*Ry.^2.*R5 );
deltaDOzz =   z1z2.*( R3 - 3*Rz.^2.*R5 );
deltaDOxy = 3*z1z2.*Rx.*Ry.*R5;
deltaDOyx = deltaDOxy;
deltaDOxz = Rx.*( z2.*R3 - 3*z1z2.*Rz.*R5 );
deltaDOzx = Rx.*( z2.*R3 + 3*z1z2.*Rz.*R5 );
deltaDOyz = Ry.*( z2.*R3 - 3*z1z2.*Rz.*R5 );
deltaDOzy = Ry.*( z2.*R3 + 3*z1z2.*Rz.*R5 );

% Rotne-Prager of reflection (factor -a^2/(4 pi eta))
DRPRxx = 1/3*R3 - Rx.^2.*R5;
DRPRyy = 1/3*R3 - Ry.^2.*R5;
DRPRzz = 1/3*R3 - Rz.^2.*R5;
DRPRxy = - Rx.Ry.*R5;
DRPRyx = DRPRxy;
DRPRxz = - Rx.Rz.*R5;
DRPRzx = DRPRxz;
DRPRyz = - Ry.Rz.*R5;
DRPRzy = DRPRyz;

% Rotne-Prager wall terms (factor a^2/(4 pi eta))
deltaDRPxx = Rz.^2.*( R5   - 5*Rx.^2.*R7 );
deltaDRPyy = Rz.^2.*( R5   - 5*Ry.^2.*R7 );
deltaDRPzz = Rz.^2.*( 3*R5 - 5*Rz.^2.*R7 );
deltaDRPxy = - 5*Rx.*Ry.*Rz.^2.*R7;
deltaDRPyx = deltaDRPxy;
deltaDRPxz = - Rx.*Rz.*( 2*R5 - 5*Rz.^2.*R7 );
deltaDRPzx = - 5*Rx.*Rz.^3.*R7;
deltaDRPyz = - Ry.*Rz.*( 2*R5 - 5*Rz.^2.*R7 );
deltaDRPzy = - 5*Ry.*Rz.^3.*R7;

%--------------------------------------------------------------------------
% Fill appropriate matrix entries via masks
%--------------------------------------------------------------------------

c4 = 3/4*diag(sigmaH);
c8 = 3/8*diag(sigmaH);
ca24 = 3/16*diag(sigmaH).^3;

D=zeros(3*nParticles,3*nParticles);

% mask to reproduce 2x2 matrix nParticles x nParticles times
dvec  = (1:3)';                          % column vector of 1:dim
dmask = dvec(:,ones(nParticles,1));     % nParticle columns each 1:dim

xxMask = logical([1 0 0; 0 0 0; 0 0 0]);
xxMask = xxMask(dmask,dmask);

yyMask = logical([0 0 0; 0 1 0; 0 0 0]);
yyMask = yyMask(dmask,dmask);

zzMask = logical([0 0 0; 0 0 0; 0 0 1]);
zzMask = zzMask(dmask,dmask);

xyMask = logical([0 1 0; 0 0 0; 0 0 0]);
xyMask = xyMask(dmask,dmask);

yxMask = logical([0 0 0; 1 0 0; 0 0 0]);
yxMask = yxMask(dmask,dmask);

xzMask = logical([0 0 1; 0 0 0; 0 0 0]);
xzMask = xzMask(dmask,dmask);

zxMask = logical([0 0 0; 0 0 0; 1 0 0]);
zxMask = zxMask(dmask,dmask);

yzMask = logical([0 0 0; 0 0 1; 0 0 0]);
yzMask = yzMask(dmask,dmask);

zyMask = logical([0 0 0; 0 0 0; 0 1 0]);
zyMask = zyMask(dmask,dmask);


D(xxMask) = -c8.*DORxx - ca24.*DRPRxx + c4.*deltaDOxx + ca24.*deltaDRPxx;

D(yyMask) = -c8.*DORyy - ca24.*DRPRyy + c4.*deltaDOyy + ca24.*deltaDRPyy;

D(zzMask) = -c8.*DORzz - ca24.*DRPRzz + c4.*deltaDOzz + ca24.*deltaDRPzz;

D(xyMask) = -c8.*DORxy - ca24.*DRPRxy + c4.*deltaDOxy + ca24.*deltaDRPxy;

D(yxMask) = -c8.*DORyx - ca24.*DRPRyx + c4.*deltaDOyx + ca24.*deltaDRPyx;

D(xzMask) = -c8.*DORxz - ca24.*DRPRxz + c4.*deltaDOxz + ca24.*deltaDRPxz;

D(zxMask) = -c8.*DORzx - ca24.*DRPRzx + c4.*deltaDOzx + ca24.*deltaDRPzx;

D(yzMask) = -c8.*DORyz - ca24.*DRPRyz + c4.*deltaDOyz + ca24.*deltaDRPyz;

D(zyMask) = -c8.*DORzy - ca24.*DRPRzy + c4.*deltaDOzy + ca24.*deltaDRPzy;
