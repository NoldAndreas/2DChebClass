function SymmetryClassificationIb_Inhomogenous()
%%  Solution I.a for a linear shear flow

%% ODE to solve
% 
% $$k(k-1)f+2y(1-k)f'+(1-y^2)f''=h(y-t)$$
%

    global dirData
    AddPaths();
    ChangeDirData([dirData filesep 'SymmetryClassification'],'ORG');

    %% Parameters
    k = 2;
    N = 100;
    L = 4;       
    
    %% Initialization
    %     
    % # for computation of similarity solution
    % # for 2D area for plotting of streamfunction
    shape     = struct('N',N,'L',L);    
    plotShape = struct('yMin',-5,'yMax',5,'N',200);

    IS         = InfSpectralLine(shape);    
    [Pts,Diff] = IS.ComputeAll(plotShape);        
    
    y    = Pts.y(2:end-1);    
    Dy   = Diff.Dy(2:end-1,2:end-1);
    DDy  = Diff.DDy(2:end-1,2:end-1);        
        
    shapeBox = struct('y1Min',-5,'y1Max',5,'N',[40,50],...
                      'y2Min',0,'y2Max',5);
    plotBox = struct('y1Min',-5,'y1Max',5,'N',[100,100],...
                     'y2Min',0,'y2Max',5);                  
                  
    BX       = Box(shapeBox);    
    [PtsBx]  = BX.ComputeAll(plotBox);
    y1       = PtsBx.y1_kv;
    y2       = PtsBx.y2_kv;
    
    zeta     = PtsBx.y1_kv./PtsBx.y2_kv;
    IP       = IS.InterpolationMatrix_Pointwise(zeta);
    
    %% Solve ODE
    %
    % # set up operator 
    % # in for-loop, set up and right hand side for time t
    % # compute and plot solution
    
    A = k*(k-1)*diag(1./(1+y.^2))+2*(1-k)*diag(y./(1+y.^2))*Dy+DDy;
    
    t = 0:0.1:3;
    f = zeros(N,length(t));
        
    for i = 1:length(t)
        h      = diag(1./(1+y.^2))*RHS(Pts.y(2:end-1)-t(i));            
        f(:,i) = [0;A\h;0];
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
    text(-4.5,-0.5,['$t= ',num2str(t(1)),'$'],'Interpreter','Latex','fontsize',18);
    text(-4.5,-0.19,['$t= ',num2str(t(ceil(end/2))),'$'],'Interpreter','Latex','fontsize',18);
    text(-4.5,-0.05,['$t= ',num2str(t(ceil(end))),'$'],'Interpreter','Latex','fontsize',18);
    
    xlabel('$x/y$','Interpreter','Latex','fontsize',25);
    ylabel('$f^{(I.b)}$','Interpreter','Latex','fontsize',25);
    
    subplot(3,2,2);
    Psi = (y2.^k).*(IP*f(:,1));        
    BX.plot(Psi,'contour',struct('clabel',false,'linecolor','k'));    
    title(['$t= ',num2str(t(1)),'$'],'Interpreter','Latex','fontsize',16);
    xlabel('$x$','Interpreter','Latex','fontsize',25);
    ylabel('$y$','Interpreter','Latex','fontsize',25);
    
    
    subplot(3,2,4);    
    Psi = (y2.^k).*(IP*f(:,ceil(end/2)));
    BX.plot(Psi,'contour',struct('clabel',false,'linecolor','k'));   
    title(['$t= ',num2str(t(ceil(end/2))),'$'],'Interpreter','Latex','fontsize',16);
    xlabel('$x$','Interpreter','Latex','fontsize',25);
    ylabel('$y$','Interpreter','Latex','fontsize',25);
    
    subplot(3,2,6);    
    Psi = (y2.^k).*(IP*f(:,end));
    BX.plot(Psi,'contour',struct('clabel',false,'linecolor','k'));   
    title(['$t= ',num2str(t(end)),'$'],'Interpreter','Latex','fontsize',16);
    xlabel('$x$','Interpreter','Latex','fontsize',25);
    ylabel('$y$','Interpreter','Latex','fontsize',25);
    
    print2eps([dirData filesep 'SelfSimilarSolution_Ib_Inh'],f1);
    saveas(f1,[dirData filesep 'SelfSimilarSolution_Ib_Inh.fig']);
    
    
    %% Homogeneous solutions
    figure('color','white');
    Psi = (y1 + y2*1i).^k;
    Psi = atan(y1./y2);
    BX.plot(Psi,'contour',struct('clabel',false,'linecolor','k'));   
    
    
    
    figure('color','white');
    for i = 1:length(t)       
        fBox  = (IP*f(:,i));
        fP    = [0;Dy*f(2:end-1,i);0];
        fPBox = (IP*fP);
        
        Psi   = (y2.^k).*fBox;
        
        title(['t = ',num2str(t(i))])
        subplot(1,2,1);
        hold off;       
        IS.plot(f(:,i));
        
        %u    = k*(y2.^(k-1)).*fBox - k*(y2.^(k-2)).*y1.*fPBox;                
        %v    = - (y2.^(k-1)).*fPBox;
        %v(y2==0) = 0;        
        %vAbs = sqrt(u.^2 + v.^2);
        subplot(1,2,2);
        BX.plot(Psi,'contour',struct('clabel',true));
        
        pause(0.1);
    end
          
    %% Right hand side of ODE
    function h = RHS(y)
        h = exp(-y.^2);
    end
end

