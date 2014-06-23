classdef ConvolutionFiniteSupport < ConvolutionFiniteSupport_NotLinear
    
    properties 
%       %There are four numberings:
%       %(1) the numberings of the domains. The points with identifyer
%       %     X are given by mark_id(:,X)       
%       mark_id
%       %(2) the numbering of the domains in the first dimension. Points
%       %in the first dimension belonging to the part with identifyer Y
%       %are given by mark_id_1(:,Y)
%       mark_id_1
%       %(3) the numbering of the domains in the second dimension. Points
%       %in the second dimension belonging to the part with identifyer Z
%       %are given by mark_id_2(:,Z)                      
%       mark_id_2
%       %mark12 makes a link between identifyers of dimensions 1,2 and 
%       % the global identifyer. Points which are in the first dimension
%       % located in a domain with identifyer Y and in the second
%       % direction in a domain with identifyer Z, belong to the global
%       % domain with the identifyer X = mark_12(Y,Z).                 
%       mark_12
    end
    
    methods (Access = public)
        X = Conv_LinearGridXY(this,ptsC,area,weights,params);       
        [X,checkSum] = Conv_LinearGridX(this,ptsC,refpts,dataAD,weights,params);        
    end
    
    methods (Abstract = true,Access = public)        
         [I1,I2] = ComputeInterpolationMatrix12(this,interp1,interp2);
         X       = ComputeConvolutionFiniteSupport(this,area,weights,pts);
    end
end