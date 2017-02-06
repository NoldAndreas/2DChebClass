classdef ComposedShape < handle
   
    properties         
        SubShape
        Int
    end
    methods
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
        function [int,area] = ComputeIntegrationVector(this)           
            int  = [];
            area = 0;
            for i = 1:length(this.SubShape)
                [intS,areaS]  = this.SubShape{i}.ComputeIntegrationVector();
                int           = [int,intS];
                if(~isempty(areaS) && ~isempty(area))
                    area          = area + areaS;
                else
                    area = [];
                end
            end            
            
            if((nargout == 1) && ~isempty(area))
                if(area == 0)
                    disp(['ComposedShape: Area is zero, Absolute error is: ',...
                                        num2str(area-sum(int))]);
                else
                    disp(['ComposedShape: Error of integration of area (ratio): ',...
                                        num2str(1-sum(int)/area)]);
                end
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
                this.SubShape{i}.ComputeInterpolationMatrix((-1:0.02:0.9)',(-1:0.02:1)',true,true);
                m = mOld + this.SubShape{i}.N1*this.SubShape{i}.N2 - 1;
                this.SubShape{i}.plot(V(mOld:m),opts);
                mOld = mOld + this.SubShape{i}.N1*this.SubShape{i}.N2;
                hold on;
            end
            hold on;            
            xlim([this.SubShape{1}.Origin(1)-5,this.SubShape{1}.Origin(1)+5]);
            xlim([this.SubShape{2}.y2Wall-5,this.SubShape{1}.Origin(2)+5]);
        end
    end    
    
end