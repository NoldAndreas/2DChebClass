classdef HalfSpaceMinusDisk < handle
   
    properties         
        SubShape
        Int
    end
    
    methods 
        function this = HalfSpaceMinusDisk(Geometry)
            %L1,R,Origin,%Rmax
            this.SubShape{1} = HalfSpaceMinusHalfDisk(Geometry);           
            
            if((Geometry.Origin(2)-Geometry.y2Wall) > Geometry.R)                
                this.SubShape{2} = StripMinusDisk(Geometry);
            elseif(((Geometry.Origin(2)-Geometry.y2Wall) < Geometry.R) ...
                            && ((Geometry.Origin(2)-Geometry.y2Wall) > 0))
                Geometry.N(2)      = ceil(Geometry.N(2)/2);
                Geometry.LeftRight = 'Left';
                this.SubShape{2}  = HalfStripMinusDisk(Geometry);
                
                Geometry.LeftRight = 'Right';
                this.SubShape{3}  = HalfStripMinusDisk(Geometry);
            elseif(((Geometry.Origin(2)-Geometry.y2Wall) < 0))            
                errror('Not yet implemented');
            end                        
        end        
        function ptsCart = GetPts(this)            
            ptsCart.y1_kv = [];
            ptsCart.y2_kv = [];            
            for i = 1:length(this.SubShape)                
                ptsCart.y1_kv = [ptsCart.y1_kv;this.SubShape{i}.Pts.y1_kv];
                ptsCart.y2_kv = [ptsCart.y2_kv;this.SubShape{i}.Pts.y2_kv];
            end            
        end                
        function ptsCart = GetCartPts(this)
            
            ptsCart.y1_kv = [];
            ptsCart.y2_kv = [];
            
            for i = 1:length(this.SubShape)
                ptsCartSub    = this.SubShape{i}.GetCartPts();    
                ptsCart.y1_kv = [ptsCart.y1_kv;ptsCartSub.y1_kv];
                ptsCart.y2_kv = [ptsCart.y2_kv;ptsCartSub.y2_kv];
            end            
        end         
        function [int] = ComputeIntegrationVector(this)           
            int = [];
            for i = 1:length(this.SubShape)
                int         = [int,this.SubShape{i}.ComputeIntegrationVector()];
            end            
        end
        
        function PlotGridLines(this)
            for i = 1:length(this.SubShape)
                this.SubShape{i}.PlotGridLines(); hold on;     
            end
        end    
        
        function PlotGrid(this)
            for i = 1:length(this.SubShape)
                this.SubShape{i}.PlotGrid(); hold on;
            end
        end    
        
        function doPlots(this,V,opts)
            if(nargin==2)
                opts = [];
            end            
            mOld = 1;
            for i = 1:length(this.SubShape)
                this.SubShape{i}.ComputeInterpolationMatrix((-1:0.02:0.6)',(-1:0.02:1)',true,true);
                m = mOld + this.SubShape{i}.N1*this.SubShape{i}.N2 - 1;
                this.SubShape{i}.doPlots(V(mOld:m),opts);
                mOld = mOld + this.SubShape{i}.N1*this.SubShape{i}.N2;
                hold on;
            end
            hold on;            
        end
    end
    
end    

