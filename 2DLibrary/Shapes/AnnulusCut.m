classdef AnnulusCut < ComposedShape
    
    methods 
        function this = AnnulusCut(Geometry)

            th             = acos(Geometry.h/Geometry.R_out);                                    
            Geometry.th1   = 3/2*pi+th;
            Geometry.th2   = 3/2*pi-th;
            this.SubShape{1}  = Wedge(Geometry);                               
            
            if(Geometry.h < Geometry.R_in)
                Geometry.leftRight = 'left';
                this.SubShape{2} = WedgeCutSide(Geometry);             
                
                Geometry.leftRight = 'right';
                this.SubShape{3} = WedgeCutSide(Geometry);             
            elseif(Geometry.h < Geometry.R_out)
                this.SubShape{2} = WedgeCut(Geometry);             
            else
                error('AnnulusCut: Define full annulus')
            end
                        
        end        
    end
end    

