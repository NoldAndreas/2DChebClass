function TimeConvergenceStudy(data)

    n1 = 5:10:25;
    n2 = 10:10:20;
    
    if(nargin == 0)    
        for i=1:length(n1)
            for j=1:length(n2)                
                data(i,j) = DiffusionAdvectionPolarInfinity(n1(i),n2(j),false);
            end 
        end
        SaveToFile('DiffusionConvergenceStudy',data,getResultsPath());
    end    
    
    times = [data.outTimes];
    error = [data.error];
    
    close all;    
    screen_size = get(0, 'ScreenSize');
    f1 = figure(1);
    set(f1, 'Position', [0 0 1000 750 ] );    
    
    set(gcf,'Color','white'); %Set background color    
    
    s       = size(times);
    %markers = {'-k',':k','-.k','om',':m','-.m'};
    markers = {'-.',':','*','.','x','o'};
    for i=1:s(2)
        plot(times(:,i),log(error(:,i)),markers{mod(i,numel(markers))+1});
                %'LineWidth',1.5);
                hold on;
    end
    legend(data.name);
    h = xlabel('t'); 
    LEG = findobj(h,'type','text');    set(LEG,'FontSize',12);
    h = ylabel('$\log_{10}(\|f_{{num}}-f_{{ana}}\|_\infty)$', 'interpreter', 'latex');
    LEG = findobj(h,'type','text');    set(LEG,'FontSize',12);

%    for i=1:length(n1)
%        for j=1:length(n2)                
%            plot(data(i,j).outTimes,log(data(i,j).error),'o');    hold on;
%        end
%    end

    


end