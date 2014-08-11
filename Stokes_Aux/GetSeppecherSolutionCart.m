function [uCart] = GetSeppecherSolutionCart(Pts,Vc,D_A,D_B,Phi)
%,p,mask_A,mask_B]
    
    if(~isstruct(Pts))
        Pts = struct('y1_kv',Pts(:,1),'y2_kv',Pts(:,2));
    end        
        
    [th,r] = cart2pol(Pts.y1_kv,Pts.y2_kv);
    
    PtsP.y1_kv = r;
    PtsP.y2_kv = th;
        
    [u]     = GetSeppecherSolution(PtsP,Vc,D_A,D_B,Phi);        %,p,mask_A,mask_B]
    [ux,uy] = GetCartesianFromPolar(u(1:end/2),u(1+end/2:end),th);
    
    uCart = [ux;uy];
    
end