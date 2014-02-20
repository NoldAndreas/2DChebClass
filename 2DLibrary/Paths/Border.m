classdef Border < Isoline
    properties
        InterpOntoBorder
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
            IP = this.InterpOntoBorder;
            
            this.IntSc      = this.IntSc*IP;
            this.IntNormal  = this.IntNormal*blkdiag(IP,IP);
            this.IntTang    = this.IntTang*blkdiag(IP,IP); 
        end     
        
        function PlotValueOnPath(this,V,ScNorTang)
            
            if(isempty(this.Interp))
               ComputeInterpolationMatrix(this,(-1:0.02:1)');
            end
            
            switch ScNorTang
                case 'normal'
                    IP = this.normal*...
                      blkdiag(this.InterpOntoBorder,this.InterpOntoBorder);                    
            end
            plot(this.Pts.t,IP*V,'o'); hold on;
            plot(this.Interp.t,this.Interp.InterPol*IP*V,'r');
        end        
        
        
    end
end