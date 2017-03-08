classdef InfCapillary_AD_FMT < InfCapillary
%OLD!!!! -> Redo and split!
   methods 
          
       function this = InfCapillary_AD_FMT(Geometry)
           this@InfCapillary(Geometry);
       end
       function [refpts,ptsy2] = GetRefY2Pts(this) 
           
           M = length(this.Pts.y1_kv);

           %Set up ptsy2
           ptsy2 = this.Pts.y2_kv;
           
           %Set up refpts
           refpts    = zeros(length(this.Pts.y1_kv),1);
           for i = 1:M
               refpts(i) = mmod(i,this.N2);               
           end                      
       end      

       function X = InterpolateAndIntegratePtsOrigin(this,ptsOr,dataIn,weights)
           
            m        = length(ptsOr.y1_kv );
            mThis    = length(this.Pts.y1_kv);
            noW      = size(weights,1);
            
            X        = zeros(m,mThis,noW+1);%always include unity weight
            checkSum = zeros(m,1);
            %Go through all points in box
            for iPts = 1:m
                
                %Get element of saved vector the distance corresponds with                

                pts = dataIn.pts;

                %center the new points aroung Pts.yi_kv(iPts)
                pts.y1_kv = pts.y1_kv + ptsOr.y1_kv(iPts);
                pts.y2_kv = pts.y2_kv + ptsOr.y2_kv(iPts);

                %Interpolate onto the new set of points
                IP    = SubShapePts(this,pts);

                %Finally get correct integration weight
                X(iPts,:,1)       = dataIn.int*IP;                
                for k = 1:noW
                    f = str2func(weights(k,:));
                    X(iPts,:,1+k) = (dataIn.int.*f(dataIn.ptsPolLoc.y2_kv)')*IP;                
                end                
                checkSum(iPts)  = dataIn.area;
            end
            
       end
              
   end
end