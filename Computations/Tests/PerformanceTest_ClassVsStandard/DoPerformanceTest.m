
function DoPerformanceTest(DC)
    %****************************************
    %************ Speed Test ****************
    % Standard structure element access
    %****************************************   
    N = 10^5;
    
    s1 = struct('ss1',struct('sss1',2));
    s2 = struct('ss2',2);
    
    c = 1;
    tic
    for i = 1:N
        c = c + s1.ss1.sss1;
    end
    toc
    
    c = 1;
    tic
    for i = 1:N
        c = c + s2.ss2;
    end
    toc    


    %****************************************
    %************ Speed Test ****************
    % Standard function evaluation
    %****************************************   
    ptsPol.y1_kv = ones(200,1)*2;
    ptsPol.y2_kv = (1:200)'/200*2*pi;
    N = 10^5;
    tic
    for i = 1:N
        ptsCart = Pol2CartPts(ptsPol);
    end
    toc
    
	tic
    for i = 1:N
        [ptsCart.y1_kv,ptsCart.y2_kv] = pol2cart(ptsPol.y2_kv,ptsPol.y1_kv);
    end
    toc
    
    
    %****************************************
    %************ Speed Test ****************
    %****************************************   
    
    if(nargin == 0)
        clear all;    
        DC  = PerformanceTestClass(2);
    end
        
    N = 10^5;
    tic
    s = struct('Ch',2*ones(N,1));
    c = 1; 
    for i = 1:N
        c = c + s.Ch(i);
    end    
    t = toc;
    disp(['Standard: ',num2str(t),' Check: ',num2str(c)]);

    tic 
    c = DC.testCount0(N);
    t = toc;
    disp(['Count in Class with struct: ',num2str(t),' Check: ',num2str(c)]);

    tic 
    c = DC.testCount(N);
    t = toc;
    disp(['Count in Class from this: ',num2str(t),' Check: ',num2str(c)]);

    tic
    c = 1;  DC.Ch = 2*ones(N,1);
    for i = 1:N
        c = c + DC.Ch(i);
    end    
    t = toc;
    disp(['Access element in obj: ',num2str(t),' Check: ',num2str(c)]);
    
    tic
    c = 1;  
    for i = 1:N
        c = c + DC.getCh(i);
    end    
    t = toc;
    disp(['Access element in obj through function: ',num2str(t),' Check: ',num2str(c)]);
    disp(['Time per function call: ',num2str(t/N)]);

    tic
    c = 1;
    Ch = 2*ones(N,1);
    for i = 1:N
        c = c+ Ch(i);
    end    
    t = toc;
    disp(['Standard, direct: ',num2str(t),' Check: ',num2str(c)]);

    tic 
    c = DC.testCount2(N);
    t = toc;
    disp(['Count in Class, direct: ',num2str(t),' Check: ',num2str(c)]);    

    %****************************************
    tic
    c = DC.sumGetCh(N);
    t = toc;
    disp(['Count in Class, calling function: ',num2str(t),' Check: ',num2str(c)]);    
    
    tic
    c = DC.sumGetChgetCh(N);
    t = toc;
    disp(['Count in Class, calling function through other function: ',num2str(t),' Check: ',num2str(c)]);    
    
    
    tic
    c = DC.sumGetChPointer(N);
    t = toc;
    disp(['Count in Class, calling function with Pointer: ',num2str(t),' Check: ',num2str(c)]);    
    
end