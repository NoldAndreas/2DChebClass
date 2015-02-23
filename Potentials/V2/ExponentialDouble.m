function [z,dzdr_r,alpha] = ExponentialDouble(r,parameter) 
%[z,dzdr_r,alpha] = Phi2DLongRange(r,parameter) 
%
%z      = -3/2*pi./((1+r.^2).^(5/2))*epsilon;
%dzdr_r = 1/r * dz/dr
%alpha  = -(pi^2/2)*epsilon = 1/2*( 2*pi*int( r*f(r), r = 0..infinity ))
    
    if(nargin == 1)
        epsilon = 1;
    elseif(isstruct(parameter))
        epsilon = parameter.epsilon;
    else
        epsilon = parameter;
    end
    lambda = 1;%parameter.lambda;

    if(isstruct(r))
        r = r.y1_kv;
    end    
    
    markG1     = (r >=1);
    markL1     = (r < 1);    
    
    z = zeros(size(r));
    
    rt         = r(markL1);
    z(markL1)  = exp(-lambda*rt.^2).*(1-erf(sqrt(-rt.^2+1)*sqrt(lambda)))*sqrt(pi/lambda);    
    
    rt         = r(markG1);
    z(markG1)  = exp(-lambda*rt.^2)*sqrt(pi/lambda);
    
    dzdr_r     = [];    
    alpha      = 4*pi*(1/(2*lambda*exp(lambda))-(1/4)*sqrt(pi)*erf(sqrt(lambda))/lambda^(3/2)+(1/4)*sqrt(pi)/lambda^(3/2));
    
    
    c      = epsilon*(-16/9*pi)/alpha;    
    z      = c*z;
    %dzdr_r = c*dzdr_r;
    alpha  = c*alpha;
    
 end