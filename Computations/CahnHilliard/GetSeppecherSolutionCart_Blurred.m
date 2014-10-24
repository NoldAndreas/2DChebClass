function [uCart] = GetSeppecherSolutionCart_Blurred(Pts,Vc,D_A,D_B,Phi)
    if(~isstruct(Pts))
        Pts = struct('y1_kv',Pts(:,1),'y2_kv',Pts(:,2));
    end        
    PtsP       = Pts;
    PtsP.y1_kv = Pts.y1_kv + 0.5;
    
    PtsM       = Pts;
    PtsM.y1_kv = Pts.y1_kv - 0.5;    

    uCart = (GetSeppecherSolutionCart(PtsM,Vc,D_A,D_B,Phi)+...
             GetSeppecherSolutionCart(Pts,Vc,D_A,D_B,Phi)+...
             GetSeppecherSolutionCart(PtsP,Vc,D_A,D_B,Phi))/3;
end