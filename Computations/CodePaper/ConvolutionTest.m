function ConvolutionTest

    set(0,'defaultaxesfontsize',20);
    set(0,'defaultlinelinewidth',2);

    defaultPos = [0 0 1000 800];

    saveDir = [pwd filesep 'Computations' filesep 'CodePaper' filesep 'Images' filesep];

    close all;

    geom.N = 50;
    geom.L = 2;  % optimize by looking at a derivative, or do parameter testing
    
    infLine = InfSpectralLine(geom);
    
    % for Gaussian plots
    PlotArea = struct('N',1000,'yMin',-20,'yMax',20);
    
    % for kernel plots
    %PlotArea = struct('N',1000,'yMin',-10,'yMax',10);
    
    [Pts,Diff,Int,Ind,Interp] = infLine.ComputeAll(PlotArea);

    % this case breaks the basic convolution
    mu1 = 0;
    mu2 = 2;   
    sigma1 = 0.25;
    sigma2 = 4;

    mu3 = 1;
    sigma3 = 1;
    
    % this case is pretty good for both (still better for ptwise
%     mu1 = 0;
%     mu2 = 1;    
%     sigma1 = 0.25;
%     sigma2 = 1;


    muConv12 = mu1+mu2;
    sigmaConv12 = sqrt(sigma1^2 + sigma2^2);

    muConv13 = mu1+mu3;
    sigmaConv13 = sqrt(sigma1^2 + sigma3^2);

    
    y = Pts.y;
    
%     yTest = y(geom.N-20);
%     shiftedKernel = Gaussian1(y-yTest);
%     figure('Position',defaultPos);
%     plot(Interp.pts,Interp.InterPol*shiftedKernel,'r');
%     hold on
%     plot(y,shiftedKernel,'or');
%     xlim([-10,10]);
%     ylim([-0.5,1.7]);
%     xlabel('$y$','interpreter','latex');
%     ylabel(['$g(y-y_0)$, $y_0=' num2str(yTest) '$'],'interpreter','latex');
%     save2pdf([saveDir 'kernel1.pdf'],gcf);
% 
%     yTest = y(geom.N-10);
%     shiftedKernel = Gaussian1(y-yTest);
%     figure('Position',defaultPos);
%     plot(Interp.pts,Interp.InterPol*shiftedKernel,'b');
%     hold on
%     plot(y,shiftedKernel,'ob');
%     xlim([-10,10]);
%     ylim([-0.5,1.7]);
%     xlabel('$y$','interpreter','latex');
%     ylabel(['$g(y-y_0)$, $y_0=' num2str(yTest) '$'],'interpreter','latex');
%     save2pdf([saveDir 'kernel2.pdf'],gcf);
% 
%     yTest = y(geom.N-6);
%     shiftedKernel = Gaussian1(y-yTest);
% 	figure('Position',defaultPos);
%     plot(Interp.pts,Interp.InterPol*shiftedKernel,'m');
%     hold on
%     plot(y,shiftedKernel,'om');
%     xlim([-10,10]);
%     ylim([-0.5,1.7]);
%     xlabel('$y$','interpreter','latex');
%     ylabel(['$g(y-y_0)$, $y_0=' num2str(yTest) '$'],'interpreter','latex');
%     save2pdf([saveDir 'kernel3.pdf'],gcf);
%     
%     return
    
    g1 = Gaussian1(y);
    g2 = Gaussian2(y);
    g3 = Gaussian3(y);
    
    gConv12 = GaussianConv12(y);
    gConv13 = GaussianConv13(y);
    
    shapeParams.N  = 200;
    shapeParams.L  = 2;
    
    % false gives pointwise convolution
    convMatrixPtwise = infLine.ComputeConvolutionMatrix(@Gaussian1,shapeParams,false);
    
    % true uses naive convolution
    convMatrixBasic = infLine.ComputeConvolutionMatrix(@Gaussian1,shapeParams,true);

    gPtwise2 = convMatrixPtwise*g2;
    gBasic2 = convMatrixBasic*g2;

    gPtwise3 = convMatrixPtwise*g3;
    gBasic3 = convMatrixBasic*g3;

    
    figure('Position',defaultPos);
    plot(Interp.pts,Interp.InterPol*g1,'r');
    hold on
    plot(Interp.pts,Interp.InterPol*g2,'b');
    plot(Interp.pts,Interp.InterPol*g3,'m');
    plot(y,g1,'or');
    plot(y,g2,'ob');
    plot(y,g3,'om');
    xlim([-20,20]);
    xlabel('$y$','interpreter','latex');
    ylabel('$g_j$','interpreter','latex');
    save2pdf([saveDir 'convolutionGaussians.pdf'],gcf);
    
    figure('Position',defaultPos);
    plot(Interp.pts,Interp.InterPol*gPtwise2,'b');
    hold on
    plot(y,gPtwise2,'ob');
    plot(Interp.pts,Interp.InterPol*gBasic2,'--b');
    plot(y,gBasic2,'xb');

    plot(Interp.pts,Interp.InterPol*gPtwise3,'m');
    hold on
    plot(y,gPtwise3,'om');
    plot(Interp.pts,Interp.InterPol*gBasic3,'--m');
    plot(y,gBasic3,'xm');

    
    %plot(y,gConv,'og');
    xlim([-20,20]);
    xlabel('$y$','interpreter','latex');
    ylabel('$g_1 \ast g_2$','interpreter','latex');
    save2pdf([saveDir 'convolutionResults.pdf'],gcf);
    
    ptwiseErr2 = L2norm(gPtwise2-gConv12)/L2norm(gConv12)
    basicErr2  = L2norm(gBasic2-gConv12)/L2norm(gConv12)
        
    ptwiseErr3 = L2norm(gPtwise3-gConv13)/L2norm(gConv13)
    basicErr3  = L2norm(gBasic3-gConv13)/L2norm(gConv13)

    
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