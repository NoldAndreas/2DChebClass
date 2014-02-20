function [X,checkSum] = Conv_NoLinearity(this,ptsC,area,weights)

    %*************************************
    %Initialization
    [convArea.int,convArea.area] = area.ComputeIntegrationVector();
    convArea.pts                 = area.GetCartPts();
    convArea.ptsPolLoc           = Cart2PolPts(convArea.pts);

    m        = length(ptsC.y1_kv);
    mThis    = length(this.Pts.y1_kv);
    noW      = length(weights);

    X        = zeros(m,mThis,noW+1);%always include unity weight
    checkSum = zeros(m,1);
    %*************************************    
    
    %Go through all points 
    for iPts = 1:m

        pts        = convArea.pts;                

        pts.y1_kv  = pts.y1_kv + ptsC.y1_kv(iPts);
        pts.y2_kv  = pts.y2_kv + ptsC.y2_kv(iPts);                                                                        

        %Interpolate onto the new set of points
        IP    = SubShapePts(this,pts);

        %Finally get correct integration weight
        X(iPts,:,1)       = convArea.int*IP;
        for k = 1:noW
            f = str2func(weights{k});
            %X(iPts,:,1+k) = (data.int.*f(data.ptsPolLoc.y2_kv)')*IP;                
            X(iPts,:,1+k) = (convArea.int.*f(convArea.ptsPolLoc)')*IP;                
            
            %f             = str2func(weights{k});
            %X(iPts,:,1+k) = ones(sum(iPts),1)*((dataDisk.int.*f(dataDisk.ptsPolLoc)')*IP);
        end
        
        
        checkSum(iPts)  = convArea.area;                                
    end
    
    
end