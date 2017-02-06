classdef Isoline < SpectralPath
    
    properties
        PhysSpace %pointer to function
        Isox1 %true if x1 is to be constant
        Backward %true if backward
        IsoValx
    end
    
    methods 
        function this = Isoline(N,ashape,isox1,backw,isovalx,polar)
            this@SpectralPath(N,polar);
            
            this.PhysSpace = @ashape.PhysSpace;
            this.Isox1     = isox1;
            this.Backward  = backw;
            this.IsoValx  = isovalx;
            
            InitializationPts(this);
        end
        
        function [y1,y2,dy1_dt,dy2_dt] = f_path(this,t)
            
            if(this.Backward)
                tt = -t;
            else
                tt = t;
            end
            O = ones(size(tt));
            
            if(this.Isox1)
                [y1,y2,J] = this.PhysSpace(this.IsoValx*O,tt);
                dy1_dt = J(:,1,2);
                dy2_dt = J(:,2,2);
            else
                [y1,y2,J] = this.PhysSpace(tt,this.IsoValx*O);
                dy1_dt = J(:,1,1);
                dy2_dt = J(:,2,1);
            end
            
            if(this.Backward)
                dy1_dt = -dy1_dt;
                dy2_dt = -dy2_dt;            
            end
            
            
        end
    end
    
end