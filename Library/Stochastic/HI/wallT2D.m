function T=wallT2D(r,sigmaH)

%--------------------------------------------------------------------------
% Extract x and z components (z normal to wall for ease of adaptation
% from 3D
%--------------------------------------------------------------------------

x=r(1:2:end);
%y=r(2:3:end);
z=r(2:2:end);

%--------------------------------------------------------------------------
% Calculate frequently-used terms
%--------------------------------------------------------------------------

nParticles=length(r)/2;

zp=z';
z1z2=z*zp;

nvec=1:nParticles;                    % row vector of column indices
nmask=nvec(ones(nParticles,1),:);     % nParticle rows of 1:nParticles (row length)

x1mx2=x(nmask)'-x(nmask);
%y1my2=y(nmask)'-y(nmask);
z1pz2=z(nmask)'+z(nmask);  % note + as we're using reflected z

% equivalent to
% z1z2=bsxfun(@times,z,z');
% x1mx2=bsxfun(@minus,x,x');
% y1my2=bsxfun(@minus,y,y');
% z1pz2=bsxfun(@plus,z,z');


% 1/|r_i-r_j|
%s1=1./( x1mx2.^2 + y1my2.^2 + z1pz2.^2).^(1/2);
s1=1./( x1mx2.^2  + z1pz2.^2).^(1/2);
s3=s1.*s1.*s1;   % multiplication faster than powers
s5=s3.*s1.*s1;


z2=zp(ones(nParticles,1),:); % z2 should be constant in columns
% equivalent to z2=repmat(z',nParticles,1); but faster

%--------------------------------------------------------------------------
% Calculate HI terms due to wall
%--------------------------------------------------------------------------

%    [-Oseen of reflection] [           wall terms          ]
Txx = - s1 - x1mx2.^2.*s3    - 2*z1z2.*(s3 - 3*s5.*x1mx2.^2);
%Tyy = - s1 - y1my2.^2.*s3    - 2*z1z2.*(s3 - 3*s5.*y1my2.^2);
Tzz = - s1 - z1pz2.^2.*s3    + 2*z1z2.*(s3 - 3*s5.*z1pz2.^2);
%Txy = - x1mx2.*y1my2.*s3     + 6*z1z2.*x1mx2.*y1my2.*s5;
%Tyx = Txy;

bracketm = z2.*s3 - 3*s5.*z1z2.*z1pz2;
bracketp = z2.*s3 + 3*s5.*z1z2.*z1pz2;

Txz = - x1mx2.*z1pz2.*s3 + 2*x1mx2.*bracketm;
Tzx = - x1mx2.*z1pz2.*s3 + 2*x1mx2.*bracketp;
%Tyz = - y1my2.*z1pz2.*s3 + 2*y1my2.*bracketm;
%Tzy = - y1my2.*z1pz2.*s3 + 2*y1my2.*bracketp;

%--------------------------------------------------------------------------
% Fill appropriate matrix entries via masks
%--------------------------------------------------------------------------

T=zeros(2*nParticles,2*nParticles);

% mask to reproduce 2x2 matrix nParticles x nParticles times
dvec=(1:2)';                          % column vector of 1:dim
dmask=dvec(:,ones(nParticles,1));     % nParticle columns each 1:dim

xxMask=logical([1 0 ; 0 0]);
xxMask=xxMask(dmask,dmask);
T(xxMask)=Txx;

% yyMask=logical([0 0 0; 0 1 0; 0 0 0]);
% yyMask=yyMask(dmask,dmask);
% T(yyMask)=Tyy;

zzMask=logical([0 0; 0 1]);
zzMask=zzMask(dmask,dmask);
T(zzMask)=Tzz;

% xyMask=logical([0 1 0; 0 0 0; 0 0 0]);
% xyMask=xyMask(dmask,dmask);
% T(xyMask)=Txy;
% 
% yxMask=logical([0 0 0; 1 0 0; 0 0 0]);
% yxMask=yxMask(dmask,dmask);
% T(yxMask)=Tyx;

xzMask=logical([0 1; 0 0]);
xzMask=xzMask(dmask,dmask);
T(xzMask)=Txz;

zxMask=logical([0 0; 1 0]);
zxMask=zxMask(dmask,dmask);
T(zxMask)=Tzx;

% yzMask=logical([0 0 0; 0 0 1; 0 0 0]);
% yzMask=yzMask(dmask,dmask);
% T(yzMask)=Tyz;
% 
% zyMask=logical([0 0 0; 0 0 0; 0 1 0]);
% zyMask=zyMask(dmask,dmask);
% T(zyMask)=Tzy;

T=3/8*diag(sigmaH)*T;