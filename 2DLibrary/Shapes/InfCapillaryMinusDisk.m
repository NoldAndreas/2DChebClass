classdef InfCapillaryMinusDisk < handle
   
    properties         
        SubShape
        Int
    end
    
    methods 
        function this = InfCapillaryMinusDisk(Geometry)
            %L1,R,Origin,%Rmax            
            
            if(Geometry.y2Max - Geometry.y2Min <= 2*Geometry.R)
                error('Not yet implemented');
            elseif(((Geometry.Origin(2)-Geometry.y2Min) >= Geometry.R) && ...
                    ((Geometry.y2Max - Geometry.Origin(2)) >= Geometry.R))
                
                Geometry.TopBottom = 'Top';
                Geometry.y2Wall    = Geometry.y2Max;
                this.SubShape{1} = StripMinusDisk(Geometry);      
                Geometry.y2Wall    = Geometry.y2Min;
                Geometry.TopBottom = 'Bottom';
                this.SubShape{2} = StripMinusDisk(Geometry);           
                
             elseif(Geometry.Origin(2) == Geometry.y2Min)                
                Geometry.TopBottom = 'Top';
                Geometry.y2Wall    = Geometry.y2Max;
                this.SubShape{1} = StripMinusDisk(Geometry);      
             elseif(((Geometry.Origin(2)-Geometry.y2Min) < Geometry.R) ...
                             && ((Geometry.Origin(2)-Geometry.y2Min) >= 0))
                         
                Geometry.TopBottom = 'Top';
                Geometry.y2Wall    = Geometry.y2Max;
                this.SubShape{1} = StripMinusDisk(Geometry);      
                        
                Geometry.TopBottom = 'Bottom';
                Geometry.N(2)      = ceil(Geometry.N(2)/2);
                Geometry.LeftRight = 'Left';
                Geometry.y2Wall    = Geometry.y2Min;
                this.SubShape{2}  = HalfStripMinusDisk(Geometry);
                 
                Geometry.LeftRight = 'Right';
                this.SubShape{3}  = HalfStripMinusDisk(Geometry);
                
            elseif(Geometry.Origin(2) == Geometry.y2Max)
                Geometry.TopBottom = 'Bottom';
                Geometry.y2Wall    = Geometry.y2Min;
                this.SubShape{1} = StripMinusDisk(Geometry);                      
                
            elseif((abs(Geometry.Origin(2)-Geometry.y2Max) < Geometry.R) ...
                            && ((Geometry.Origin(2)-Geometry.y2Max) <= 0))
                
                Geometry.TopBottom = 'Bottom';
                Geometry.y2Wall    = Geometry.y2Min;
                this.SubShape{1} = StripMinusDisk(Geometry);      
                         
                Geometry.TopBottom = 'Top';
                Geometry.N(2)      = ceil(Geometry.N(2)/2);
                Geometry.LeftRight = 'Left';
                Geometry.y2Wall    = Geometry.y2Max;
                this.SubShape{2}  = HalfStripMinusDisk(Geometry);
                 
                Geometry.LeftRight = 'Right';
                this.SubShape{3}  = HalfStripMinusDisk(Geometry);
            else
                error('Not yet implemented');
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
        
        function plot(this,V,opts)
            if(nargin==2)
                opts = [];
            end            
            mOld = 1;
            for i = 1:length(this.SubShape)
                this.SubShape{i}.ComputeInterpolationMatrix((-1:0.02:0.6)',(-1:0.02:1)',true,true);
                m = mOld + this.SubShape{i}.N1*this.SubShape{i}.N2 - 1;
                this.SubShape{i}.plot(V(mOld:m),opts);
                mOld = mOld + this.SubShape{i}.N1*this.SubShape{i}.N2;
                hold on;
            end
            hold on;            
        end
    end
    
end    

