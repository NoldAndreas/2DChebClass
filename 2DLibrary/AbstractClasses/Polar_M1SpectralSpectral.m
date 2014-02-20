classdef Polar_M1SpectralSpectral < M1SpectralSpectral        
    
    methods
        
        function this = Polar_M1SpectralSpectral(N1,N2)
             this@M1SpectralSpectral(N1,N2);                          
             this.polar = 'polar';
        end              
        
        function Int = ComputeIntegrationVector(this,t1Odd)
            if(nargin > 1)
                ComputeIntegrationVector@M1SpectralSpectral(this,t1Odd);
            else
                ComputeIntegrationVector@M1SpectralSpectral(this);
            end
            Int             = this.Int.*this.Pts.y1_kv'; 
            Int(Int == inf) = 0; 
            Int(isnan(Int)) = 0;
            this.Int = Int;
        end        
        
        function Diff = ComputeDifferentiationMatrix(this)
            Diff  = ComputeDifferentiationMatrix@M1SpectralSpectral(this);
            Diff  = PolarDiffOperators(Diff,this.Pts);           
        end              

    end
    
end