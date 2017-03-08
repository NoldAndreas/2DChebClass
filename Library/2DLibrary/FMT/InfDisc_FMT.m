classdef InfDisc_FMT < InfDisc & ConvolutionFiniteSupport_NotLinear
    methods (Access = public)
        
         function this = InfDisc_FMT(Geometry)   
            
             this@InfDisc(Geometry);
            
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
                         
        function AD = ComputeConvolutionFiniteSupport(this,area,weights,pts)
            %only tested if pts are in a shape of a HalfSpace
            %INPUT:
            % pts - structure with 'y1_kv','y2_kv','y1','y2'
            
            fprintf('Computing interpolation for matrices for averaged densities..\n');
            tic            
            AD = Conv_NoLinearity(this,pts,area,weights);                
            disp([num2str(toc),'s']);  
        end        
        function [AD,AAD] = GetAverageDensities(this,area,weights)
            %**********************************************************
            %AD  - Average densities to get average densities
            %AAD - average the average densities to compute free energy
            %
            %**********************************************************
            AD   = ComputeConvolutionFiniteSupport(this,area,weights,this.Pts);
            AAD  = AD;
        end   

        
    end
    
end