classdef Wedge_Interface < Polar_M1SpectralSpectral

    properties   
        R_in = 0;
        R_out
        th1,th2                         
        sphere
    end
    
    methods        
        function this = Wedge_Interface(Geometry)
            this@Polar_M1SpectralSpectral(Geometry.N(1),Geometry.N(2));
            
            if(isfield(Geometry,'y1Min'))
                this.R_in  = Geometry.y1Min;
                this.R_out = Geometry.y1Max;
                this.th1   = Geometry.y2Min;
                this.th2   = Geometry.y2Max;
                
            else
                
                if(isfield(Geometry,'R'))
                    this.R_out    = Geometry.R;
                else
                    this.R_out    = Geometry.R_out;
                end
                if(isfield(Geometry,'R_in'))
                        this.R_in     = Geometry.R_in;
                end
                this.th1      = Geometry.th1;
                this.th2      = Geometry.th2;  
            end
            
            InitializationPts(this);  
            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;                
            end            
            if(isfield(Geometry,'sphere'))
                this.sphere = Geometry.sphere;
            end
        end   
    end
    
    methods (Access = public)

        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)
            [y1_kv,y2_kv,J,dH1,dH2] = M1TrefPolarWall(x1,x2,this.R_out,pi/2,0.8);                                    
        end
        
        function [x1,x2] = CompSpace(this,y1,y2)             
             error = MException('InfCapillary_Interface',...
                                'inverse function not yet implemented');
             throw(error);
        end             
        
    end
end