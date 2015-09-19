function [AAD] = ComputeConvolutionFiniteSupport(this,area,weights,pts)
%GetAverageDensities(this,area,weights,pts)            
    %AAD - average the average densities to compute free energy
    fprintf('Computing interpolation for matrices for averaging of averaged densities..\n');    
    
    y2Lim         = this.Sub_Strip.y2Max + 0.5/sin(this.alpha);
    markSub       = (pts.y2_kv < y2Lim);
    ptsSub        = pts;
    ptsSub.y1_kv  = pts.y1_kv(markSub);
    ptsSub.y2_kv  = pts.y2_kv(markSub);
    ptsSub.x1_kv  = pts.x1_kv(markSub);
    ptsSub.x2_kv  = pts.x2_kv(markSub);
    ptsSub.y2     = pts.y2(pts.y2 < y2Lim);
    
    
    AAD_Substrip              = zeros(length(pts.y1_kv),this.Sub_Strip.M,length(weights)+1);
    X                         = this.Sub_Strip.ComputeConvolutionFiniteSupport(area,weights,ptsSub);
    AAD_Substrip(markSub,:,:) = X;
    AAD_SubHS                 = this.Sub_HalfSpace.ComputeConvolutionFiniteSupport(area,weights,pts);
    
    AAD = zeros(size(AAD_SubHS,1),...
                size(AAD_Substrip,2)+size(AAD_SubHS,2),...
                size(AAD_SubHS,3));
    for i = 1:size(AAD_SubHS,3)
        AAD(:,:,i) = [AAD_Substrip(:,:,i),AAD_SubHS(:,:,i)];
    end
    
    %AAD  = Conv_LinearGridXY(this,pts,area,weights);
end