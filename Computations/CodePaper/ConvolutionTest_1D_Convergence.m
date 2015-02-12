function ConvolutionTest_1D_Convergence

    close all;

    geom.N = 30;
    geom.L = 2;  % optimize by looking at a derivative, or do parameter testing
	PlotArea = struct('N',1000,'yMin',-20,'yMax',20);
    
    infLine = InfSpectralLine(geom);
    [Pts,Diff,Int,Ind,Interp] = infLine.ComputeAll(PlotArea);

    % this case breaks the basic convolution
    mu1 = 0;
    mu2 = 2;   
    sigma1 = 0.25;
    sigma2 = 4;

    mu3 = 1;
    sigma3 = 1;

    muConv12 = mu1+mu2;
    sigmaConv12 = sqrt(sigma1^2 + sigma2^2);

    muConv13 = mu1+mu3;
    sigmaConv13 = sqrt(sigma1^2 + sigma3^2);

    
    y              = Pts.y;           
    shapeParams.L  = 2;
    %shapeParamsHalfLim.L    = 2;
    shapeParamsHalfLim.yMin = 1;
    shapeParamsHalfLim.yMax = 10;
    shapeParamsLim.yMin = -1;
    shapeParamsLim.yMax = 1;
    
    % false gives pointwise convolution
    NS_d = 5;
    NS = 20:NS_d:80;%20:NS_d:40;
    for i = 1:length(NS)
        shapeParams.N     = NS(i);
        shapeParamsHalfLim.N  = NS(i);
        shapeParamsLim.N  = NS(i);
        res{i}.ConvBH1_HalfLim  = infLine.ComputeConvolutionMatrix(@BH1,shapeParamsHalfLim,false);
        res{i}.ConvGauss1       = infLine.ComputeConvolutionMatrix(@Gaussian1,shapeParams,false);
        %res{i}.ConvBH1     = infLine.ComputeConvolutionMatrix(@BH1,shapeParams,false);        
        res{i}.ConvBH1_Lim      = infLine.ComputeConvolutionMatrix(@BH1,shapeParamsLim,false);
        res{i}.NS         = NS(i);
    end        
    
    cols = {'b','m','k','r','g'}; nosyms = 5;
    syms = {'d','s','o','>','<'}; nocols = 5;
    legendstring = {};
    
	figure('color','white','Position',[0 0 800 800]);         
    PlotMatrixErrorOverY('ConvGauss1','o','k','Gaussian');
    %PlotMatrixErrorOverY('ConvBH1','o','m','BH1');
    PlotMatrixErrorOverY('ConvBH1_HalfLim','o','g','BH1_{HalfLim}');
    PlotMatrixErrorOverY('ConvBH1_Lim','o','b','BH1_{Lim}');

    xlabel('$N$','Interpreter','Latex','fontsize',15);        
    set(gca,'YScale','log');
    set(gca,'linewidth',1.5);
    set(gca,'fontsize',15);
    legend(legendstring,'Location','southoutside','Orientation','horizontal');        
    
     function PlotMatrixErrorOverY(A_name,sym,col,name)     
        leg_string = {};        
        
        y = Pts.y;
        n = 1;
        for j = 1:(length(res)-1)            
            line(n)   = max(max(abs(res{j}.(A_name)-res{j+1}.(A_name))));
            line_N(n) = (res{j}.NS);%+res(k1,k2+1).NS)/2;                
            n = n+1;                
        end       
        plot(line_N,line,['-',sym,col],'MarkerSize',10,'MarkerFaceColor',col); hold on;
        legendstring(end+1) = {name};
    end
    
    
    % true uses naive convolution
    convMatrixBasic = infLine.ComputeConvolutionMatrix(@Gaussian1,shapeParams,true);

    gPtwise2 = convMatrixPtwise*g2;
    gBasic2 = convMatrixBasic*g2;

    gPtwise3 = convMatrixPtwise*g3;
    gBasic3 = convMatrixBasic*g3;

     
%     figure('Position',defaultPos);
%     plot(Interp.pts,Interp.InterPol*g1,'r');
%     hold on
%     plot(Interp.pts,Interp.InterPol*g2,'b');
%     plot(Interp.pts,Interp.InterPol*g3,'m');
%     plot(y,g1,'or');
%     plot(y,g2,'ob');
%     plot(y,g3,'om');
%     xlim([-20,20]);
%     xlabel('$y$','interpreter','latex');
%     ylabel('$g_j$','interpreter','latex');
%     save2pdf([saveDir 'convolutionGaussians.pdf'],gcf);
%     
%     figure('Position',defaultPos);
%     plot(Interp.pts,Interp.InterPol*gPtwise2,'b');
%     hold on
%     plot(y,gPtwise2,'ob');
%     plot(Interp.pts,Interp.InterPol*gBasic2,'--b');
%     plot(y,gBasic2,'xb');
% 
%     plot(Interp.pts,Interp.InterPol*gPtwise3,'m');
%     hold on
%     plot(y,gPtwise3,'om');
%     plot(Interp.pts,Interp.InterPol*gBasic3,'--m');
%     plot(y,gBasic3,'xm');
% 
%     
%     %plot(y,gConv,'og');
%     xlim([-20,20]);
%     xlabel('$y$','interpreter','latex');
%     ylabel('$g_1 \ast g_2$','interpreter','latex');
%     save2pdf([saveDir 'convolutionResults.pdf'],gcf);
    
    ptwiseErr2 = L2norm(gPtwise2-gConv12)/L2norm(gConv12)
    basicErr2  = L2norm(gBasic2-gConv12)/L2norm(gConv12)
        
    ptwiseErr3 = L2norm(gPtwise3-gConv13)/L2norm(gConv13)
    basicErr3  = L2norm(gBasic3-gConv13)/L2norm(gConv13)

    function g = BH1(x)
        x = abs(x);
        %g = x.*BarkerHenderson_2D(x); g(x==inf) = 0; g(x==-inf) = 0;
        %g(x <= 1) = BarkerHenderson_2D(1);
        g = 1./(x.^4);
        %g(x <= 1) = BarkerHenderson_2D(1);
        g(x <= 1) = 1;
    end
    
    function g = Gaussian1(x)
        g = 1/sqrt(2*pi)/sigma1 * exp( -(x-mu1).^2/(2*sigma1^2) );
    end

    function g = Gaussian2(x)
        g = 1/sqrt(2*pi)/sigma2 * exp( -(x-mu2).^2/(2*sigma2^2) );
    end

    function g = Gaussian3(x)
        g = 1/sqrt(2*pi)/sigma3 * exp( -(x-mu3).^2/(2*sigma3^2) );
    end


    function g = GaussianConv12(x)
        g = 1/sqrt(2*pi)/sigmaConv12 * exp( -(x-muConv12).^2/(2*sigmaConv12^2) );
    end

    function g = GaussianConv13(x)
        g = 1/sqrt(2*pi)/sigmaConv13 * exp( -(x-muConv13).^2/(2*sigmaConv13^2) );
    end


    function Norm = L2norm(f)
        Norm = Int*(f.^2);
    end


end