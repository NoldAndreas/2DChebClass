classdef InfCapillary < InfCapillaryGeneral & ConvolutionFiniteSupport

	methods        
        function this = InfCapillary(Geometry)
            this@InfCapillaryGeneral(Geometry);
            
            %There are four numberings: 
             % (here, there is only one domain,simplifying the numbering, as
             % all identifiers = 1)
             %(1) the numberings of the domains. The points with identifyer
             %     X are given by mark_id(:,X)       
             this.mark_id = true(this.M,1);
             %(2) the numbering of the domains in the first dimension. Points
             %in the first dimension belonging to the part with identifyer Y
             %are given by mark_id_1(:,Y)
             this.mark_id_1 = true(this.N1,1);
             %(3) the numbering of the domains in the second dimension. Points
             %in the second dimension belonging to the part with identifyer Z
             %are given by mark_id_2(:,Z)                      
             this.mark_id_2 = true(this.N2,1);
             %mark12 makes a link between identifyers of dimensions 1,2 and 
             % the global identifyer. Points which are in the first dimension
             % located in a domain with identifyer Y and in the second
             % direction in a domain with identifyer Z, belong to the global
             % domain with the identifyer X = mark_12(Y,Z).                 
             this.mark_12   = ones(1,1); 
        end
    end
    
    methods (Access = public)
        function [y1,dy1,dx,ddx,dddx,ddddx] = PhysSpace1(this,x1)
            [y1,dy1,dx,ddx,dddx,ddddx] = SqrtMap(x1,this.L1,inf);            
            y1 =  this.y10 + y1; 
        end
        function xf = CompSpace1(this,y1)            
            xf  = InvSqrtMap(y1 - this.y10,this.L1,inf);
        end
        
        AD = ComputeConvolutionFiniteSupport(this,area,weights,pts,params)                 
    end
end