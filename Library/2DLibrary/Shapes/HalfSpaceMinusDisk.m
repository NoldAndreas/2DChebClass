classdef HalfSpaceMinusDisk < ComposedShape
   
    methods 
        function this = HalfSpaceMinusDisk(Geometry)
            %L1,R,Origin,%Rmax
            this.SubShape{1} = HalfSpaceMinusHalfDisk(Geometry);           
            
            if((Geometry.Origin(2)-Geometry.y2Wall) > Geometry.R)                
                this.SubShape{2} = StripMinusDisk(Geometry);
            elseif(((Geometry.Origin(2)-Geometry.y2Wall) < Geometry.R) ...
                            && ((Geometry.Origin(2)-Geometry.y2Wall) > 0))
                Geometry.N(2)      = ceil(Geometry.N(2)/2);
                Geometry.LeftRight = 'Left';
                Geometry.TopBottom = 'Bottom';
                this.SubShape{2}  = HalfStripMinusDisk(Geometry);
                
                Geometry.LeftRight = 'Right';
                this.SubShape{3}  = HalfStripMinusDisk(Geometry);
            elseif(((Geometry.Origin(2)-Geometry.y2Wall) < 0))            
                errror('Not yet implemented');
             end                        
        end  
    end
end    

