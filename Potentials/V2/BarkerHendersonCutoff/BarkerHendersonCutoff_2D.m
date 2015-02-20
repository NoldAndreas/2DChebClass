function [z,dzdr_r,alpha] = BarkerHendersonCutoff_2D(r,parameter) 
%[z,dzdr_r,alpha] = Phi2DLongRange(r,parameter) 
%
%z      = f(r);
%dzdr_r = 0
%alpha  =/2*( 2*pi*int( r*f(r), r = 0..infinity ))

    global r_cutoff    
    
    if(nargin < 2)
        epsilon = 1;
    elseif(isstruct(parameter))
        epsilon = parameter.epsilon;
        r_cutoff = parameter.r_cutoff;
    else
        epsilon = parameter;
    end
        
    if(isstruct(r))
        r = r.y1_kv;
    end
        
    markG1 = (r <= r_cutoff) & (r >= 1);      
    markL1 = (r < 1);
    
    rG1    = r(markG1);
    rL1    = r(markL1);
    
    z            = 0*r;    
    z(markG1)    = 2*BH_2D_I(rG1,sqrt(r_cutoff^2-rG1.^2));
    z(markL1)    = BH_2D_L1(rL1);
    dzdr_r       = 0*r;
    alpha        = 1/2*( BH_0D_I(r_cutoff) - BH_0D_I(1) );
    
    %Rescale
    c = (-16/9*pi)/alpha;    
    z      = c*z;
    dzdr_r = c*dzdr_r;
    alpha  = c*alpha;

    %Rescale with epsilon
    z      = epsilon*z;
    alpha  = epsilon*alpha;
end