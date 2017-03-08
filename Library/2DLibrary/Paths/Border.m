classdef Border < Isoline
    properties
        InterpOntoBorder
        
        IntSC_Path,IntNormal_Path,IntTang_Path
    end
    methods        
        function this = Border(N,ashape,border,polar)        
            switch border
                case 'left'
                    isox1 = true;
                    backw = true;
                    isovalx = ashape.Pts.x1(1);
                case 'right'
                    isox1 = true;
                    backw = false;
                    isovalx = ashape.Pts.x1(end);
                case 'top'
                    isox1 = false;
                    backw = true;
                    isovalx = ashape.Pts.x2(end);
                case 'bottom'
                    isox1 = false;
                    backw = false;
                    isovalx = ashape.Pts.x2(1);
            end
            this@Isoline(N,ashape,isox1,backw,isovalx,polar);
            this.InterpOntoBorder = ashape.SubShapePts(this.Pts);
        end                
        
        function ComputeIntegrationVector(this)            
            ComputeIntegrationVector@SpectralPath(this);
            
            this.IntSC_Path      = this.IntSc; 
            this.IntNormal_Path  = this.IntNormal;
            this.IntTang_Path    = this.IntTang_Path;
            
            IP              = this.InterpOntoBorder;
            
            this.IntSc      = this.IntSc*IP;
            this.IntNormal  = this.IntNormal*blkdiag(IP,IP);
            this.IntTang    = this.IntTang*blkdiag(IP,IP); 
        end     
        
        function PlotValueOnPath(this,V,ScNorTang,var)
            
            if(isempty(this.Interp))
               ComputeInterpolationMatrix(this,(-1:0.02:1)');
            end
            
            if((nargin >= 3) && strcmp(ScNorTang,'normal'))
                IP = this.normal*...
                      blkdiag(this.InterpOntoBorder,this.InterpOntoBorder);
            else
                IP = this.InterpOntoBorder;
            end
            
            if((nargin < 4) || (nargin >= 4) && strcmp(var,'t'))
                pt_t     = this.Pts.t;
                interp_t = this.Interp.t;
            elseif(strcmp(var,'y1'))
                pt_t     = this.Pts.y1_kv;
                interp_t = this.Interp.pts_y1;
            elseif(strcmp(var,'y2'))
                pt_t     = this.Pts.y2_kv;
                interp_t = this.Interp.pts_y2;
            end
            plot(pt_t,IP*V,'o'); hold on;
            plot(interp_t,this.Interp.InterPol*IP*V,'r');
        end        
        
        
    end
end