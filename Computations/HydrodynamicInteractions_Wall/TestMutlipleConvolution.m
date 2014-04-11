function TestMutlipleConvolution
%%
% 
% $$I = f({\bf r}) \int_{|{\bf r} - \tilde {\bf r}| > 1} f(\tilde{\bf r})
% W({\bf r} - \tilde {\bf r}) d\tilde {\bf r}$$
%
% with
%
% $$W(r) = \log r$$
% 
% For convenience this is computed as $I = I_1 - I_2$ where
% 
% $$I_1 = f({\bf r}) \int f(\tilde{\bf r}) \bar W({\bf r} - \tilde {\bf r}) d\tilde {\bf r}$$
%
% and
%
% $$I_{2} = f({\bf r}) \int_{|{\bf r} - \tilde {\bf r}| < 1} f(\tilde{\bf r}) \bar W({\bf r} - \tilde {\bf r}) d\tilde {\bf r}$$
% 
% where we extend W(r) to be valid for $r<1$, by
%
% $$ \bar W(r) = W(r)\left( 1 - e^{- \left( \frac{r}{1-r} \right)^2  } \right) $$
%
    

    PhysArea = struct('N',[20,20],...
                      'L1_Skewed',4,...
                      'L2_Skewed',2,...
                      'L2_AD_Skewed',2.,...
                      'y2wall',0.,...
                      'N2bound',14,'h',1,...
                      'alpha_deg',90);
                  
	PlotArea = struct('y1Min',-5,'y1Max',15,'y2Min',1,'y2Max',15,'N1',80,'N2',80);

    HS                 = HalfSpaceSkewed_FMT(PhysArea,1);
	[Pts,Diff,Int,Ind] = HS.ComputeAll();
	PtsCart            = HS.GetCartPts();
    HS.InterpolationPlotCart(PlotArea,true);
    
    M = length(PtsCart.y1_kv);
    
    %% Initialize Convolution
    ConvMod = zeros(M,M);
    for i = 1:M
        ConvMod(i,:) = Int.*(HIFilm_2D(PtsCart.y1_kv(i),PtsCart.y2_kv(i))).';
    end            
    
    shapeDisc.R   = 1;
    shapeDisc.N   = [20,20];
    area          = Disc(shapeDisc);
    ConvDisc      = HS.ComputeConvolutionFiniteSupport(area,{'HI_Weight_Test'},HS.Pts);
    
    areaPlot = struct('y1Min',0,'y1Max',1,'y2Min',0,'y2Max',2*pi,'N1',80,'N2',80);
    area.ComputeAll(areaPlot);
    
    figure('Name','Weight on Disc to be subtracted');
    area.doPlots(HI_Weight_Test(area.Pts));       
    title('Weight on Disc to be subtracted');
    
    ft            = f(PtsCart);
    
    %% Test 1 - 11/04/2014    
    
    ConvFull = ConvMod - ConvDisc(:,:,2);
    figure('Name','Result of full convolution I');
    HS.doPlots(ft.*(ConvFull*ft));              
    
    
    %% Auxiliary function 
    
    function y = f(pts)
        y = exp(- ((pts.y1_kv).^2 + (pts.y2_kv).^2)/20);
    end

    function y = HIFilm_2D(y1,y2)
        
        pts.y1_kv      = sqrt((PtsCart.y1_kv-y1).^2 + (PtsCart.y2_kv-y2).^2);     
        pts.y1_kv(isnan(pts.y1_kv)) = 0;        
        y = HI_Weight_Test(pts);        
        
    end
    
   % Test 2    
    % 	shapeAnn.RMin = 1;
%     shapeAnn.L    = 1;
%     shapeAnn.N    = [20,20];
%     area          = InfAnnulus(shapeAnn);    
   % ConvAnn       = HS.ComputeConvolutionFiniteSupport(area,{'RotnePragner_f1','RotnePragner_f2','RotnePragner_f0'},HS.Pts);

%     subplot(2,2,1);    
%     HS.doPlots(ft);
%     subplot(2,2,2);    
%     HS.doPlots(ConvAnn(:,:,2)*ft);
%     subplot(2,2,3);    
%     HS.doPlots(ConvAnn(:,:,3)*ft);
%     subplot(2,2,4);    
%     HS.doPlots(ConvAnn(:,:,4)*ft);


end