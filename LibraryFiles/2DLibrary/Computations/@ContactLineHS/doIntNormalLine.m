function I = doIntNormalLine(this,y2Max,y1,f_loc,f_hs)
    
    
    if(~isempty(this.IntMatrFex))
        %Integrate in three intervals: 
        %(1) [0 0.5]
        [pts.y2_kv,w] = getPointsWeight(20,[0 0.5]);
        pts.y1_kv     = y1*ones(size(pts.y2_kv));
        IP            = this.HS.AD.SubShapePtsCart(pts);    

        I = w*IP*f_hs;

        %(2) [0.5 1]
        [pts.y2_kv,w] = getPointsWeight(20,[0.5 1]);
        pts.y1_kv     = y1*ones(size(pts.y2_kv));

        IP_hs         = this.HS.AD.SubShapePtsCart(pts);
        IP_loc        = this.HS.SubShapePtsCart(pts);

        I   = I + w*(IP_hs*f_hs + IP_loc*f_loc);

        %(3) [1 y2Max]
        [pts.y2_kv,w] = getPointsWeight(400,[1 y2Max]);
        pts.y1_kv     = y1*ones(size(pts.y2_kv));

        IP_hs         = this.HS.AD.SubShapePtsCart(pts);
        IP_loc        = this.HS.SubShapePtsCart(pts);

        I   = I + w*(IP_hs*f_hs + IP_loc*f_loc);        
    else
        %Integrate in one intervals: 
        [pts.y2_kv,w] = getPointsWeight(500,[0.5 y2Max]);
        pts.y1_kv     = y1*ones(size(pts.y2_kv));        
        IP            = this.HS.SubShapePtsCart(pts);
        I             = w*IP*(f_hs + f_loc);        
    end
        
    function [yp,w] = getPointsWeight(n,y2I)        
        [xCheb,wCheb] = ClenCurtFlip(n);
        w             = wCheb*(y2I(2)-y2I(1))/2;    
        yp            = y2I(1) + (1+xCheb)*(y2I(2)-y2I(1))/2;
    end
end