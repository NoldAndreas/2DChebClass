function SymmetryClassificationIb_Inhomogenous2()
%%  Solution I.a for a linear shear flow

%% ODE to solve
% 
% $$k(k-1)f+2y(1-k)f'+(1-y^2)f''=h(y-t)$$
%

    global dirData
    AddPaths();
    ChangeDirData([dirData filesep 'SymmetryClassification'],'ORG');

    %% Parameters
    k = 1;
    N = 200;         
    t = 0:0.05:0.1;
    
    %% Initialization
    %     
    % # for computation of similarity solution
    % # for 2D area for plotting of streamfunction
    shape     = struct('N',N,'yMin',-5,'yMax',5);
    plotShape = struct('yMin',-5,'yMax',5,'N',200);

    IS         = SpectralLine(shape);    
    [Pts,Diff] = IS.ComputeAll(plotShape);        
    
    y    = Pts.y;%(2:end-1);    
    Dy   = Diff.Dy;%(2:end-1,2:end-1);
    DDy  = Diff.DDy;%(2:end-1,2:end-1);        
        
    shapeBox = struct('y1Min',-1,'y1Max',1,'N',[41,50],...
                      'y2Min',0,'y2Max',1);
    plotBox = struct('y1Min',-1,'y1Max',1,'N',[100,100],...
                     'y2Min',0,'y2Max',1);  
                  
    BX       = Box(shapeBox);    
    [PtsBx]  = BX.ComputeAll(plotBox);
    y1       = PtsBx.y1_kv;
    y2       = PtsBx.y2_kv;
    
    zeta     = PtsBx.y1_kv./PtsBx.y2_kv;
    zeta(zeta>max(y)) = max(y);
    zeta(zeta<min(y)) = min(y);
    zeta(PtsBx.y1_kv == 0) = 0;
    IP       = IS.InterpolationMatrix_Pointwise(zeta);    
    
    
    Geometry = struct('R_in',0.,'R_out',1,...
                      'th1',0,'th2',pi,'N',[10,20]);
    
    WD        = Wedge(Geometry);
    WD.ComputeAll();
    WD.ComputeInterpolationMatrix((-1:0.02:1)',(-1:0.02:1)',true,true);
    PtsWD = WD.GetCartPts();
    xWD   = PtsWD.y1_kv;  yWD = PtsWD.y2_kv;    
    IPWD = BX.SubShapePtsCart(PtsWD);
    
    %% Solve ODE
    %
    % # set up operator 
    % # in for-loop, set up and right hand side for time t
    % # compute and plot solution
    
    A = k*(k-1)*diag(1./(1+y.^2))+2*(1-k)*diag(y./(1+y.^2))*Dy+DDy;        
    f = zeros(N,length(t));
        
    IntM = zeros(N);
	vh   = zeros(1,N);
    for i = 2:N
        %Matrix for integral int(v(y),y=y..Y)
        hh          = y(i) - y(i-1);
        vh([i-1,i]) = vh([i-1,i]) + hh/2;
        IntM(i,:)    = vh;
    end
    
    for i = 1:length(t)
        
        g1 = RHS(y-t(i))./(1+y.^2);        
        f(:,i) = y.*(IntM*g1)-IntM*(g1.*y);% +sqrt(pi)/4;        
        f(:,i) = f(:,i) - (f(1,i)+f(end,i))/2;
        
%         h      = diag(1./(1+y.^2))*RHS(y-t(i));                    
%         %h      = diag(1./(1+Pts.y.^2))*RHS(Pts.y-t(i));                    
%         %h      = RHS(Pts.y(2:end-1)-t(i)); 
%         f0   = sqrt(pi)/4;
%         fend = -sqrt(pi)/4;
%         f(:,i) = [f0;
%                   DDy(2:end-1,2:end-1)\...
%                     (h(2:end-1)-(f0*A(2:end-1,1) + fend*A(2:end-1,end)));...
%                  fend];
%        %f(:,i) = Diff.DDy\h;
    end
    
    %% Plot    
    % 
	% # Plot
    % # Video
    % 
    
    f1 = figure('color','white','Position',[0 0 1000 900]);
    
    subplot(3,2,[1,3,5]);
    IS.plot(f(:,1),'plain'); hold on;
    IS.plot(f(:,ceil(end/2)),'plain');
    IS.plot(f(:,ceil(end)),'plain');
    text(-4.5,0.23,['$t= ',num2str(t(1)),'$'],'Interpreter','Latex','fontsize',18);
    text(-4.5,0.3,['$t= ',num2str(t(ceil(end/2))),'$'],'Interpreter','Latex','fontsize',18);
    text(-4.5,0.35,['$t= ',num2str(t(ceil(end))),'$'],'Interpreter','Latex','fontsize',18);
    
    xlabel('$x/y$','Interpreter','Latex','fontsize',25);
    ylabel('$f^{(I.b)}$','Interpreter','Latex','fontsize',25);
    
    subplot(3,2,2);    
    Psi = (y2.^k).*(IP*f(:,1));        
    Psi(y2.^k == inf) = 0;       
    BX.plot(Psi,'contour',struct('clabel',false,'linecolor','k'));    hold on;
    u = BX.Diff.Dy2*Psi;    v = -BX.Diff.Dy1*Psi;
    WD.plotFlux([IPWD*u;IPWD*v],[],[],1.5,'k');
    %plotFlux(this,flux,maskAdd,fl_norm,lw,c,plain)
    xlim([-1 1]); ylim([0 1]);
    title(['$t= ',num2str(t(1)),'$'],'Interpreter','Latex','fontsize',20);
    xlabel('$x$','Interpreter','Latex','fontsize',20);
    ylabel('$y$','Interpreter','Latex','fontsize',20);
    
    
    subplot(3,2,4);    
    Psi = (y2.^k).*(IP*f(:,ceil(end/2)));
    Psi(y2.^k == inf) = 0;
    BX.plot(Psi,'contour',struct('clabel',false,'linecolor','k'));      hold on; 
    u = BX.Diff.Dy2*Psi;    v = -BX.Diff.Dy1*Psi;
    WD.plotFlux([IPWD*u;IPWD*v],[],[],1.5,'k');
    xlim([-1 1]); ylim([0 1]);
    title(['$t= ',num2str(t(ceil(end/2))),'$'],'Interpreter','Latex','fontsize',20);
    xlabel('$x$','Interpreter','Latex','fontsize',20);
    ylabel('$y$','Interpreter','Latex','fontsize',20);
    
    subplot(3,2,6);    
    Psi = (y2.^k).*(IP*f(:,end));
    Psi(y2.^k == inf) = 0;
    BX.plot(Psi,'contour',struct('clabel',false,'linecolor','k'));   hold on;
    u = BX.Diff.Dy2*Psi;    v = -BX.Diff.Dy1*Psi;
    WD.plotFlux([IPWD*u;IPWD*v],[],[],1.5,'k');
    xlim([-1 1]); ylim([0 1]);
    title(['$t= ',num2str(t(end)),'$'],'Interpreter','Latex','fontsize',20);
    xlabel('$x$','Interpreter','Latex','fontsize',20);
    ylabel('$y$','Interpreter','Latex','fontsize',20);
    
    print2eps([dirData filesep 'SelfSimilarSolution_Ib_Inh'],f1);
    saveas(f1,[dirData filesep 'SelfSimilarSolution_Ib_Inh.fig']);
    
    
%     %% Homogeneous solutions
%     figure('color','white');
%     %Psi = (y1 + y2*1i).^k;
%     Psi = atan(y1./y2);
%     BX.plot(Psi,'contour',struct('clabel',false,'linecolor','k'));   
    
    figure('color','white');
    for i = 1:length(t)       
        fBox  = (IP*f(:,i));
        fP    = Dy*f(:,i);
        %[0;Dy*f(2:end-1,i);0];
        fPBox = (IP*fP);
        
        Psi   = (y2.^k).*fBox;
        Psi(y2.^k == inf) = 0;
        
        title(['t = ',num2str(t(i))])
        subplot(1,2,1);
        hold off;       
        IS.plot(f(:,i));
        
        %u    = k*(y2.^(k-1)).*fBox - k*(y2.^(k-2)).*y1.*fPBox;                
        %v    = - (y2.^(k-1)).*fPBox;
        %v(y2==0) = 0;        
        %vAbs = sqrt(u.^2 + v.^2);
        subplot(1,2,2);
        BX.plot(Psi,'contour',struct('clabel',true)); hold on;
        %BX.plot(Psi);%,'contour',struct('clabel',true));
        %figure; 
        BX.plotFlux([BX.Diff.Dy2*Psi,-BX.Diff.Dy1*Psi]);
        
        pause(0.1);
    end
          
    %% Right hand side of ODE
    function h = RHS(y)
        %h = (1+y.^2).*y.*exp(-y.^2);
        h = y.*exp(-y.^2);
    end
end

