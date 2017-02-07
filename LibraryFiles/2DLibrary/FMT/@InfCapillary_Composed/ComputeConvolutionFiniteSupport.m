function [AAD] = ComputeConvolutionFiniteSupport(this,area,weights,pts)
%GetAverageDensities(this,area,weights,pts)            
    %AAD - average the average densities to compute free energy
    %%fprintf('Computing interpolation for matrices for averaging of averaged densities..');    
    %%AAD  = Conv_LinearGridXY(this,pts,area,weights);        
    
    %GetAverageDensities(this,area,weights,pts)            
    %AAD - average the average densities to compute free energy
    fprintf('Computing interpolation for matrices for averaging of averaged densities..\n');    
    
    alpha            = pi/2;
    
    markBottom       = (pts.y2_kv < this.Bottom_Strip.y2Max + 0.5/sin(alpha));
    ptsBottom        = pts;
    ptsBottom.y1_kv  = pts.y1_kv(markBottom);
    ptsBottom.y2_kv  = pts.y2_kv(markBottom);
    ptsBottom.x1_kv  = pts.x1_kv(markBottom);
    ptsBottom.x2_kv  = pts.x2_kv(markBottom);
    ptsBottom.y2     = pts.y2(pts.y2 < this.Bottom_Strip.y2Max + 0.5/sin(alpha));    
    
    
    markTop       = (pts.y2_kv > this.Top_Strip.y2Min - 0.5/sin(alpha));
    ptsTop        = pts;
    ptsTop.y1_kv  = pts.y1_kv(markTop);
    ptsTop.y2_kv  = pts.y2_kv(markTop);
    ptsTop.x1_kv  = pts.x1_kv(markTop);
    ptsTop.x2_kv  = pts.x2_kv(markTop);
    ptsTop.y2     = pts.y2(pts.y2 > this.Top_Strip.y2Min - 0.5/sin(alpha));
    
    AAD_BottomStrip                = zeros(length(pts.y1_kv),this.Bottom_Strip.M,length(weights)+1);
    AAD_TopStrip                   = zeros(length(pts.y1_kv),this.Top_Strip.M,length(weights)+1);
    
    AAD_BottomStrip(markBottom,:,:) = this.Bottom_Strip.ComputeConvolutionFiniteSupport(area,weights,ptsBottom);
    AAD_MainStrip                  = this.Main_Strip.ComputeConvolutionFiniteSupport(area,weights,pts);
    AAD_TopStrip(markTop,:,:)      = this.Top_Strip.ComputeConvolutionFiniteSupport(area,weights,ptsTop);        
    
    AAD = zeros(size(AAD_BottomStrip,1),...
                size(AAD_BottomStrip,2)+size(AAD_MainStrip,2)+size(AAD_TopStrip,2),...
                size(AAD_BottomStrip,3));
    for i = 1:size(AAD_BottomStrip,3)
        AAD(:,:,i) = [AAD_BottomStrip(:,:,i),AAD_MainStrip(:,:,i),AAD_TopStrip(:,:,i)];
    end
    
    %AAD  = Conv_LinearGridXY(this,pts,area,weights);
end  