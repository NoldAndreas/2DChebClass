function [X,checkSum] = InterpolateAndIntegratePtsOrigin(this,ptsOr,data,weights)

    m        = length(ptsOr.y1_kv);
    mThis    = length(this.Pts.y1_kv);
    noW      = size(weights,1);

    X        = zeros(m,mThis,noW+1);
    checkSum = zeros(m,1);

    %Go through all points in box
    for iPts = 1:m

        j   = mod(iPts,ptsOr.N2);
        if(j==0) 
            j = ptsOr.N2; 
        end 
        %*******************************
        % Origin:    in Box
        % Integrate: in Box
        %Get element of saved vector the distance corresponds with                
        if(data(j).area ~= -1)
            pts = data(j).pts;

            %center the new points aroung Pts.yi_kv(iPts)
            pts.y1_kv = pts.y1_kv + ptsOr.y1_kv(iPts);

            %Interpolate onto the new set of points
            IP    = SubShapePts(this,pts);

            %Finally get correct integration weight
            X(iPts,:,1)       = data(j).int*IP;                
            for k = 1:noW
                f = str2func(weights(k,:));
                X(iPts,:,1+k) = (data(j).int.*f(data(j).ptsPolLoc.y2_kv)')*IP;                
            end                                        
            checkSum(iPts)  = data(j).area;
        end
    end
end