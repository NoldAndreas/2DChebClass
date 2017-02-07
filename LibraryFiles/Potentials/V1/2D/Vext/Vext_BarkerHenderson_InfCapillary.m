function [VBack_S,VAdd_S]=Vext_BarkerHenderson_InfCapillary(y1,y2,t,opts)
%Currently designed for HalfSpace

% stoc: inputs are (x,[],t,opts)
% ddft: inputs are (r,theta,t,opts)

    y2Min       = opts.y2Min;
    y2Max       = opts.y2Max;
    epw         = opts.epsilon_w;            

    if(~isfield(opts,'epsilon_w_end'))
        epw_end = epw;
    else
        epw_end = opts.epsilon_w;
    end

    if(isfield(opts,'tau'))
        t       = t/opts.tau;
    end

	VBack       = zeros(size(y1));        
    DVBackDy1   = zeros(size(y1));
    DVBackDy2   = zeros(size(y1));

    VBack_S = struct('V',VBack,...
                'dy1',DVBackDy1,'dy2',DVBackDy2,...
                'grad',[DVBackDy1;DVBackDy2]);            
    
    [wall_1,dwall_1] = BarkerHendersonWall(y2 - y2Min,1);
    [wall_2,dwall_2] = BarkerHendersonWall(y2Max - y2,1);

    VAdd    = wall_1*(epw(1) + (epw_end(1) - epw(1))*(1 -exp(-t^2))) ...
                + wall_2*(epw(2) + (epw_end(2) - epw(2))*(1 -exp(-t^2)));
    dVAdd   = dwall_1*(epw(1) + (epw_end(1) - epw(1))*(1 -exp(-t^2))) ...
                + dwall_2*(epw(2) + (epw_end(2) - epw(2))*(1 -exp(-t^2)));
    
    VAdd_S      = struct('V',VAdd);
    VAdd_S.dy1  = zeros(size(VAdd));
    VAdd_S.dy2  = dVAdd;
    
end