function RecordDDFT(input)
    frames = 10;
    %*****************************
    %Initialization of Data     
    v2struct(input);
    v2struct(data);
    
    [N1,N2,h1,h2,y1Plot,y2Plot,plotTimes] = LoadNumData(optsNum);    
    
    if(isfield(data,'Subspace'))
        v2struct(Subspace);
    end
    
    xl = [(min(y1Plot)-0.5) (max(y1Plot)+0.5)];
    yl = [(min(y2Plot)-0.5) (max(y2Plot)+0.5)];

    
    ma3 = (~Ind.bound  & (Pts.y1_kv >=  min(y1Plot)) & (Pts.y1_kv <=  max(y1Plot)) & (Pts.y2_kv <=  max(y2Plot)) & (Pts.y2_kv >=  min(y2Plot)));
    maB = (Ind.bound   & (Pts.y1_kv >=  min(y1Plot)) & (Pts.y1_kv <=  max(y1Plot)) & (Pts.y2_kv <=  max(y2Plot)) & (Pts.y2_kv >=  min(y2Plot)));
    maC = (Ind.corners & (Pts.y1_kv >=  min(y1Plot)) & (Pts.y1_kv <=  max(y1Plot)) & (Pts.y2_kv <=  max(y2Plot)) & (Pts.y2_kv >=  min(y2Plot)));       
    
    fl_norm = 0.1*max(max(max(abs(flux_t)))); xs = min(xl); ys = min(yl);
    %**************************************
    %Initialization of figure and screen, and movie
    close all;
    figure
    
    gifFile = getMovieFile('Movie');            
    f1 = figure(1);
    screen_size = get(0, 'ScreenSize');
    %set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
    set(f1, 'Position', [0 0 550 850]);    
    set(gcf,'Color','white'); %Set background color    
    
       
    d      = ceil(length(plotTimes)/frames);
    
    %*****************************
    k = 1; fileNames = [];
    for i=1:d:length(plotTimes)        
        x       = X_t(:,i);
        rho     = rho_t(:,i);
        t       = plotTimes(i);
                
%         Plot Flux-Vectorgraph
            hold off;
            fl    = flux_t(:,i);
            fl_y1 = fl(1:N1*N2);
            fl_y2 = fl(N1*N2+1:end);
            y1_s  = Pts.y1_kv;  y2_s  = Pts.y2_kv;               
        
            NormQuiverPlot(y1_s(ma3),y2_s(ma3),fl_y1(ma3),fl_y2(ma3),fl_norm,xs,ys); hold on;                
            doPlots_IP_Contour(data.Interp,rho);        hold on;         
            xlim(xl); ylim(yl);
            xlabel('y1');        ylabel('y2');
            title(['t = ',num2str(t)]);
            pbaspect([(xl(2)-xl(1)) (yl(2)-yl(1)) 1]);
        
            %Record(i,gifFile);getPDFMovieFile('Movie',i)
            fileName = getPDFMovieFile('Movie1',k);
            %fileName = ['Movie1',num2str(k),'.pdf'];
            save2pdf(fileName,gcf);
            k = k+1;
            fileNames = [fileNames,' ',fileName];
    end       
    
    str = [getMovieFile([],true),optsNum.DDFTCode,'_' getTimeStr()]; %returns minute as char
    allPdfFiles = [str,'.pdf'];
    swfFile     = [str,'.swf'];
    
    system(['C:\pdftk.exe ', fileNames ,' cat output ',allPdfFiles]);    
    system(['C:\pdf2swf.exe -s framerate=7 -o ',swfFile,' ', allPdfFiles]); 
    system(['copy ',getPDFMovieFile('Movie1',1),' ',str,'POSTER.pdf']);
    system(['del ',fileNames]);
end