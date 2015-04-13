function [VBack_S,VAdd_S]=Vext_SquareWell(y1,y2,t,opts)
%Currently designed for HalfSpace

% stoc: inputs are (x,[],t,opts)
% ddft: inputs are (r,theta,t,opts)

    lambdaW         = 0.5;
    epsilon_w       = opts.epsilon_w;        
    
    VBack       = zeros(size(y1));
    DVBackDy1   = zeros(size(y1));
    DVBackDy2   = zeros(size(y1));

    VBack_S = struct('V',VBack,...
                'dy1',DVBackDy1,'dy2',DVBackDy2,...
                'grad',[DVBackDy1;DVBackDy2]);                    
            
	VAdd        = -epsilon_w*((y2-0.5)<=lambdaW);
    VAdd_S      = struct('V',VAdd);
    VAdd_S.dy1  = zeros(size(y1));
    VAdd_S.dy2  = [];%zeros(size(y1));

end