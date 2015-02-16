classdef AnnulusCut < handle
   
    properties         
        WD,WC                
    end
    
    methods 
        function this = AnnulusCut(Geometry)

            th             = acos(Geometry.h/Geometry.R_out);                                    
            Geometry.th1   = 3/2*pi+th;
            Geometry.th2   = 3/2*pi-th;
            this.WD        = Wedge(Geometry);                               
            
            this.WC        = WedgeCut(Geometry);             
                        
        end        
        function ptsCart = GetCartPts(this)   
            ptsCart = struct('y1_kv',[this.WD.GetCartPts.y1_kv;...
                                      this.WC.GetCartPts.y1_kv],...
                             'y2_kv',[this.WD.GetCartPts.y2_kv;...
                                     this.WC.GetCartPts.y2_kv]);                                                
        end         
        function [int,area] = ComputeIntegrationVector(this)           
            [WDint,WD_a] = this.WD.ComputeIntegrationVector();                        
            [WCint,WC_a] = this.WC.ComputeIntegrationVector();
            
            int       = [WDint,WCint];
            area      = WD_a + WC_a;  
            
            if(area == 0)
                disp(['AnnulusCut: Area is zero, Absolute error is: ',...
                                    num2str(area-sum(int))]);
            else
                disp(['AnnulusCut: Error of integration of area (ratio): ',...
                                    num2str(1-sum(int)/area)]);
            end
        end
        function PlotGridLines(this)
            this.WD.PlotGridLines(); hold on;
            this.WC.PlotGridLines();
        end
        function PlotGrid(this)
            this.WD.PlotGrid(); hold on;
            this.WC.PlotGrid();
        end
%         function plot(this,V)
%             this.WD.plot(V(1:this.WD.M));
%             this.WC.plot(V(1+this.WD.M:end));
%         end
    end
    
end    

