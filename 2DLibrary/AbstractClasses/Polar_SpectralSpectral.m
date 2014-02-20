classdef Polar_SpectralSpectral < SpectralSpectral
    
    methods
        
        function this = Polar_SpectralSpectral(N1,N2)
            this@SpectralSpectral(N1,N2);
        end
        function Int = ComputeIntegrationVector(this,t1Odd,t2Odd)
            if(nargin == 1)
                ComputeIntegrationVector@SpectralSpectral(this);
            elseif(nargin == 2)
                ComputeIntegrationVector@SpectralSpectral(this,t1Odd);
            elseif(nargin == 3)
                ComputeIntegrationVector@SpectralSpectral(this,t1Odd,t2Odd);
            end
            Int             = this.Int.*this.Pts.y1_kv'; 
            Int(Int == inf) = 0; 
            Int(isnan(Int)) = 0;
            this.Int = Int;
            this.polar       = 'polar';
        end        
        
        function Diff = ComputeDifferentiationMatrix(this)
            Diff  = ComputeDifferentiationMatrix@SpectralSpectral(this);
            Diff  = PolarDiffOperators(Diff,this.Pts);           
            this.Diff = Diff;
        end              

    end
    
end