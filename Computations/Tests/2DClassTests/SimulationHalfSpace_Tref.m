function data = SimulationHalfSpace_Tref()

    disp('** Simulation HalfSpace Tref **');
    AddPaths();    
    close all;
    
    %Initialization
    if(nargin == 0)
        PhysArea = struct('L1',2,'L2',2,'N',[20,20]);
        PlotArea = struct('y1Min',-6,'y1Max',6,'N1',100,...
                          'y2Min',0,'y2Max',10,'N2',70);                
        vext  = @Vext17;
    end
    
    THS                       = TrefHalfSpace(PhysArea);
    [Pts,Diff,Int,Ind,Interp] = THS.ComputeAll(PlotArea);
            
    [V,Vdiff]          = vext(Pts.y1_kv,Pts.y2_kv);    
    [VP,VPDiff]        = vext(Interp.pts1,Interp.pts2);  
    
    
    THS.doPlots(V,'SC');    
    figure
    subplot(2,2,1); THS.doPlots(Vdiff.dy1,'SC');
    subplot(2,2,2); THS.doPlots(Vdiff.dy2,'SC');
    subplot(2,2,3); THS.doPlots(Vdiff.ddy1,'SC');
    subplot(2,2,4); THS.doPlots(Vdiff.ddy2,'SC');
    
	figure;
    subplot(2,2,1); THS.doPlots(Diff.Dy1*V,'SC');
    subplot(2,2,2); THS.doPlots(Diff.Dy2*V,'SC');
    subplot(2,2,3); THS.doPlots(Diff.DDy1*V,'SC');
    subplot(2,2,4); THS.doPlots(Diff.DDy2*V,'SC');   
    
    %Check Differentiation
    vplot    = Interp.InterPol*V;        
    dataOrg  = displayErrorsPos(Pts,vplot,VP,V,Vdiff,Diff,'cart');
    
    %***********************************************************
    %Update Poles
    %***********************************************************
    close all;
    disp('Updated Poles');
    [Vt,Pts,Diff,Int,Ind,Interp] = THS.UpdatePadeValues(V,PlotArea);
                                            
    [V,Vdiff]   = vext(Pts.y1_kv,Pts.y2_kv);    
    [VP]        = vext(Interp.pts1,Interp.pts2);                          
    
    disp(['Error of interpolation:',num2str(max(abs(V-Vt)))]);
        
    figure;
    THS.doPlots(V,'SC'); 
    THS.PlotLineOfPoles(V);
    
    displayErrorsPos(Pts,Interp.InterPol*V,VP,V,Vdiff,Diff,'cart');
    
    %***********************************************************
    %Update Poles
    %***********************************************************
    disp('Updated Poles');
    [Vt,Pts,Diff,Int,Ind,Interp] = THS.UpdatePadeValues(V,PlotArea);
    
    [V,Vdiff]   = vext(Pts.y1_kv,Pts.y2_kv);    
    [VP]        = vext(Interp.pts1,Interp.pts2);                          
    disp(['Error of interpolation:',num2str(max(abs(V-Vt)))]);
    
    figure;
    THS.doPlots(V,'SC'); 
    THS.PlotLineOfPoles(V);    

    
    disp('Updated Poles');
    [Vt,Pts,Diff,Int,Ind,Interp] = THS.UpdatePadeValues(V,PlotArea);

    [V,Vdiff]   = vext(Pts.y1_kv,Pts.y2_kv);    
    [VP]        = vext(Interp.pts1,Interp.pts2);              
    disp(['Error of interpolation:',num2str(max(abs(V-Vt)))]);
    
    figure;
    THS.doPlots(V,'SC'); 
    THS.PlotLineOfPoles(V);    

    
    disp('Updated Poles');
    [Vt,Pts,Diff,Int,Ind,Interp] = THS.UpdatePadeValues(V,PlotArea);

    
    [V,Vdiff]   = vext(Pts.y1_kv,Pts.y2_kv);    
    [VP]        = vext(Interp.pts1,Interp.pts2);                          
    disp(['Error of interpolation:',num2str(max(abs(V-Vt)))]);    
    
    figure;
    THS.doPlots(V,'SC'); 
    THS.PlotLineOfPoles(V);    
    
    figure;
    subplot(2,2,1); THS.doPlots(Vdiff.dy1,'SC');
    subplot(2,2,2); THS.doPlots(Vdiff.dy2,'SC');
    subplot(2,2,3); THS.doPlots(Vdiff.ddy1,'SC');
    subplot(2,2,4); THS.doPlots(Vdiff.ddy2,'SC');
    
    figure;
    subplot(2,2,1); THS.doPlots(Diff.Dy1*V,'SC');
    subplot(2,2,2); THS.doPlots(Diff.Dy2*V,'SC');
    subplot(2,2,3); THS.doPlots(Diff.DDy1*V,'SC');
    subplot(2,2,4); THS.doPlots(Diff.DDy2*V,'SC');    
     
    %Check Differentiation 
    data    = displayErrorsPos(Pts,Interp.InterPol*V,VP,V,Vdiff,Diff,'cart');
           
    %Check Integration
%    data.Int = abs(Int*V-VInt);
%    display([' Error in Integration: ', num2str(data.Int)]);                
        
    %Check Convolution
    %fP = f1(Pts.y1_kv,Pts.y2_kv);
    %data.Conv = max(abs(Interp.InterPol*(Conv*fP) - fConv(Interp.pts1,Interp.pts2)));
    %display([' Error in Convolution: ', num2str(data.Conv)]);
    
%     data.N1 = N1; data.N2 = N2;
% %    data.Interp = Interp; data.f = V;    
%     
%     %******** Plotting **********
%     figure
%     set(gcf,'Color','white'); %Set background color    
%     
%     subplot(1,2,1);
%     TS.doPlots(V);    
%     title('Interpolation');    
%     pbaspect([1 1 1]);
%     
%     subplot(1,2,2);
%     TS.doPlots(fConv(Interp.pts1,Interp.pts2));
%     title('Convolution');
%     pbaspect([1 1 1]);    
    
end