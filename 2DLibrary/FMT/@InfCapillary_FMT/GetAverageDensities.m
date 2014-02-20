function [AD,AAD] = GetAverageDensities(this,area,weights)
%**********************************************************
%AD  - Average densities to get average densities
%AAD - average the average densities to compute free energy
%
%**********************************************************
    [dataS.int,dataS.area] = area.ComputeIntegrationVector();
    dataS.pts              = area.GetCartPts();
    dataS.ptsPolLoc        = Cart2PolPts(dataS.pts);

    [refpts,ptsy2]         = this.AD.GetRefY2Pts(this.y2wall + 2*this.R);

    if(isa(area,'Disc') && area.sphere)
        for iPts = 1:length(ptsy2)
            dataAD(iPts) = Intersect_Disk(this,ptsy2(iPts),area.R,[area.N1,area.N2],'sphere');
        end
    elseif(isa(area,'Disc') && ~area.sphere)
        for iPts = 1:length(ptsy2)
            dataAD(iPts) = Intersect_Disk(this,ptsy2(iPts),area.R,[area.N1,area.N2]);
        end             
     elseif(isa(area,'Ball'))
        for iPts = 1:length(ptsy2)
            dataAD(iPts) = Intersect_Ball(this,ptsy2(iPts),area.R,[area.N1,area.N2]);
        end
    else
        exc = MException('HalfSpace_FMT:GetAverageDensities','case not implemented');
        throw(exc);
    end

    fprintf('Computing interpolation for matrices for averaged densities..');
    tic
    [AD,checkAD] = InterpolateAndIntegratePtsOrigin(this,refpts,dataS,dataAD,weights);            
    t = toc;
    disp([num2str(t),'s']);           
    
	AAD = this.AD.GetAverageDensities(area,weights,this.Pts);
    %**********************************************************            
    %**********************************************************            
    %Test:
    [errAD,ierrAD] = max(abs(checkAD - sum(AD(:,:,1),2)));
    y1err = this.AD.Pts.y1_kv(ierrAD);
    y2err = this.AD.Pts.y2_kv(ierrAD);            

    disp(['Max. Error in AD: ', num2str(errAD),...
                  ' at y_1= ', num2str(y1err),...
                  ' at y_2= ', num2str(y2err)]);                      
    %**********************************************************                       

end   