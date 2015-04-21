function [x,w] = ClenCurtFlip(N)
% [x,w] = ClenCurt(N) takes in the number of discretisation points N and outputs x
% and w, being the usual Chebyshev points in [1,-1] and w, the Clenshaw-Curtis
% integration weights (same as Trefethan weights in clencurt.m)

% Compute Chebychev collocation points
if(N==0)
    x = 1;
    w = 1;
    return;
end

theta = pi*(0:N)'/N;
x = flipdim(cos(theta),1);

w    = zeros(1,N+1);    % Initialise w, the weights
ii   = 2:N;             
v    = ones(N-1,1);     

% For even numbers
if mod(N,2)==0          % Check for even total cheb points
    w(1)    = 1/(N^2-1);   
    w(N+1)  = w(1);     
    for k=1:N/2-1       
        v = v - 2*cos(2*k*theta(ii))/(4*k^2 - 1);   
    end
    v = v - cos(N*theta(ii))/(N^2-1);               
else
    w(1)   = 1/N^2;      
    w(N+1) = w(1);     
    for k=1:(N-1)/2    
        v = v - 2*cos(2*k*theta(ii)) / (4*k^2 - 1);
    end
end

w(ii) = 2*v/N;        

