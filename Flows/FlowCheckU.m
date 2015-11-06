function U_S=FlowCheckU(y1S,y2S,t,optsPhys)

U0 = optsPhys.U0;

%--------------------------------------------------------------------------   

%theta = mod(atan(y2S./y1S),2*pi)
%r = sqrt(y1S.^2 + y2S.^2);

[theta,r] = cart2pol(y1S,y2S);

U1        = -r.*sin(theta);
U2        = r.*cos(theta);
DU1Dy1    = zeros(size(y1S)); %NOT CORRECT
DU2Dy2    = zeros(size(y1S)); %NOT CORRECT

%U1        = U0*ones(size(y1S));
%U2        = zeros(size(y1S));
%DU1Dy1    = zeros(size(y1S));
%DU2Dy2    = zeros(size(y1S));

U_S = struct('U',[U1;U2],'U1',U1,'U2',U2,'div',DU1Dy1+DU2Dy2);

end