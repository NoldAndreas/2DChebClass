function TrefClosedBox2D_BurgersEquation()
%************************************************
%
% Solves 
%   (DYN 1) du/dt = Lap(u)+e^u
%   (BC)    u     = 0
%************************************************
    global dirData
    disp(' ** Adaptive Grid Burgers Equation **');
    close all;
    
    %************************************************
    %****************  Preprocess  ****************
    %************************************************        
    nu          = 0.01; %0.001
    
    N1          = 15; 
    N2          = 40;
    PlotArea    = struct('y1Min',-1,'y1Max',1,'N1',50,...
                         'y2Min',-1,'y2Max',1,'N2',50);
              
    TB                         = BoxTrefSpectralSpectral(N1,N2);
    [Pts,Diff,Int,Ind,Interp]  = TB.ComputeAll(PlotArea);
    
    %***********************************************
    %*************** Initial Condition *************
    %***********************************************       
    Dt           = 0.01;
    t_n          = 0;
        
    x1h          = (1+Pts.y1_kv)/2;    
    x2h          = (1+Pts.y2_kv)/2;    
    A            = 0.5+x1h;
    u_n          = A.*(sin(2*pi*x2h) + 0.5*sin(pi*x2h));
                
	TB.plot(u_n,'SC');
    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
	mM            = ones(N1*N2,1);
    mM(Ind.bound) = 0;
    opts          = odeset('RelTol',10^-8,'AbsTol',10^-8,'Mass',diag(mM));   
    
    
    fileNames = [];
    
    k = 1;
    for k = 1:54
        fileName = getPDFMovieFile('Movie1',k);
       fileNames = [fileNames,' ',fileName];
    end
    
    while(t_n < 0.3)
        [u_n,t_n] = odeSolver(u_n,t_n,Dt);
        %[u_n,t_n] = EulerForward(u_n,t_n,Dt);
        hold off;
        
        TB.plot(u_n,'SC');
        title(['t = ',num2str(t_n)]);
        
        %Adapt Grid
        if(t_n>=0.1)
            %Dt = 0.2;
            [u_n,Pts,Diff,Int,Ind,Interp] = TB.UpdatePadeValues(u_n,PlotArea);
            TB.PlotLineOfPoles(u_n);
            Dt  = 0.1/max(abs(ft(t_n,u_n)));
            disp(['Dt = ',num2str(Dt)]);          
        end        
        view([1 1 1]);
        zlim([-2 2]);
        %pause(0.05);
        
        fileName = getPDFMovieFile('Movie1',k);
        %fileName = ['Movie1',num2str(k),'.pdf'];
        save2pdf_new(fileName,gcf);
        k = k+1;
        fileNames = [fileNames,' ',fileName];           
         
        %close all;
        
    end          
    
    % Save Movie
    str         = [dirData filesep 'Version1'];
    allPdfFiles = [str,'.pdf'];
    swfFile     = [str,'.swf'];

    %system(['C:\pdftk.exe ', fileNames ,' cat output ',allPdfFiles]);    
    switch computer
		case {'MAC','MACI','MACI64'}			
            system(['/usr/local/bin/pdftk ', fileNames ,' cat output ',allPdfFiles]);    
		case {'PCWIN','PCWIN64'}
            system(['C:\pdftk ', fileNames ,' cat output ',allPdfFiles]);    
        otherwise
            gs= 'gs';
    end
    
    system(['C:\pdf2swf.exe -s framerate=5 -o ',swfFile,' ', allPdfFiles]);
    system(['copy ',getPDFMovieFile('Movie1',1),' ',str,'POSTER.pdf']);
    system(['del ',fileNames]);       
    disp(['Swf Movie` saved in: ',swfFile]);
    
    %*********************************************
    function dudt = ft(t,u)               
        
        dudt            = nu*Diff.Lap*u-(diag(u)*Diff.Dy2)*u;

        dudt(Ind.bound) = u(Ind.bound);
        dudt(Ind.right) = Ind.normalRight*Diff.grad*u;
        dudt(Ind.left)  = Ind.normalLeft*Diff.grad*u;        
    end

    function [u_nP1,t_nP1]=odeSolver(u_nM1,t_nM1,Dt)                
        [outTimes,u_t]   = ode15s(@ft,[t_nM1,t_nM1+Dt],u_nM1,opts);
        t_nP1            = outTimes(end);
        u_nP1            = u_t(end,:)';
    end

end


% 
%         % for swf Recording
%         fileName = getPDFMovieFile('Movie1',k);
%         %fileName = ['Movie1',num2str(k),'.pdf'];
%         save2pdf(fileName,gcf);
%         k = k+1;
%         fileNames = [fileNames,' ',fileName];           
%         
%         close all;
%     end
% 
%     %% Save Movie
%     str         = [dirData filesep 'Equilibrium' filesep this.FilenameEq ,'_DensitySlices'];
%     allPdfFiles = [str,'.pdf'];
%     swfFile     = [str,'.swf'];
% 
%     system(['C:\pdftk.exe ', fileNames ,' cat output ',allPdfFiles]);    
%     system(['C:\pdf2swf.exe -s framerate=5 -o ',swfFile,' ', allPdfFiles]);
%     system(['copy ',getPDFMovieFile('Movie1',1),' ',str,'POSTER.pdf']);
%     system(['del ',fileNames]);       
%     disp(['Swf Movie` saved in: ',swfFile]);