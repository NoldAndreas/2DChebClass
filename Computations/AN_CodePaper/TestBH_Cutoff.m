function TestBH_Cutoff

    r_cutoff = 2.5;

    f = @BarkerHendersonCutoff_2D; %BarkerHendersonHardCutoff_2D
    
    optsPhys.r_cutoff = r_cutoff;
    optsPhys.epsilon  = 1;

    shape.yMin = 1;
    shape.yMax = r_cutoff;
    shape.N    = 20;    
    plotRange  = struct('yMin',1,'yMax',r_cutoff,'N',100);
    
    Int   = SpectralLine(shape);
    Int.ComputeAll(plotRange);
    
    %Check interpolation
    fx = Phi(Int.Pts.y);
    
    
    PrintErrorPos(max(abs(Int.Interp.InterPol - Phi(Int.InterPol.pts))),' Interpolation error');
    Interval.Interp.pts
    
    function z = Phi(x)
       z = f(x,optsPhys); 
    end
    
end