function [z,DiffInt] = M1TrefInterval(x,y1,y2,y0,ep_y)

    L = (y2-y1)/2;
    
    d  = InvLinearMap(y0,y1,y2);
    ep = ep_y/L;        
    
    [y,Diff] = M1Tref(x,d,ep);        
    z        = LinearMap(y,y1,y2);
    
    DiffInt.dydx   = Diff.dydx*L;
    DiffInt.dyddx  = Diff.dyddx*L;
    
    DiffInt.dydy0  = Diff.dyd_d;
    DiffInt.dyddy0 = Diff.dydd_d;
    
    DiffInt.dydxdy0 = Diff.dydxd_d*L;
    
    DiffInt.dydep_y = Diff.dyde;
    DiffInt.dyddep_y = Diff.dydde;
    DiffInt.dydxdep_y = Diff.dydedx*L;

end