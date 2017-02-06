function [y,dy] = HaleTee(x)
    
z = 0.1; n = 3;

    m = 0.523231225073770;
    z2k = [0.830135290736502 + 0.557562013657515i
    0.221599693267731 + 0.975137721526374i];


    m14 = m^(.25); % 4th root of elliptic paramter
    L = -2*log(m14)/pi;
    
    %original(for evaluation of complex arguments): ellipkkp,ellipjc
    % from Driscoll’s Schwarz–Christoffel toolbox
    % this implementation is for real arguments only
    
    [K Kp] = ellipke(L); % elliptic integrals 
    h1 = m14*ellipj(2*K*asin(z(:))/pi,L); % the map to the disk (2.3)
    % evaluate the summation on the rhs of (2.18)
    ZZ = repmat([-m14;m14;h1],1,n-1); Z2K = repmat(z2k.',length(h1)+2,1);
    ZZ1 = ZZ - Z2K; idx1 = find(real(ZZ1)<0 & imag(ZZ1)>=0);
    WW1 = log(ZZ1); WW1(idx1) = WW1(idx1) - 2i*pi;
    ZZ2 = ZZ - conj(Z2K); idx2 = find(real(ZZ2)<0 & imag(ZZ2)<0);
    WW2 = log(ZZ2); WW2(idx2) = WW2(idx2) + 2i*pi;
    
    ak = diff(d)/pi; % ak given by jumps in delta
    sumlogs = 1i*((WW1 - WW2)*ak); % the summation
    % system of equations for A, a0, b0
    M = (1-m14)^2/(4*m14); M = [2/(m14^2-1) .5 .5 ; -1 -M 1+M ; 1 -1-M M];
    rhs = [ak'*(angle(z2k)-pi)+d(1) ; [-1;1]-sumlogs(1:2)];
    lhs = -(1-m14^2)/(1+m14^2)*(M*rhs);
    A = lhs(1); a0 = lhs(2); b0 = lhs(3);
    g = A + a0./(h1-1) + b0./(h1+1) + sumlogs(3:end); % the map g(z) = h3(h1(z))
end