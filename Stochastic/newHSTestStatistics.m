function newHSTestStatistics(dt,q0,v,D,N)


    function p = p1(q,dt,q0,v,D)
        p = 1/sqrt(4*pi*D*dt) * exp( - (q-q0-v*dt).^2/(4*D*dt) );
    end
     
    function p = p2(q,dt,q0,v,D)
        p = exp(-v*q0/D) / sqrt(4*pi*D*dt) .* exp( - (q + q0 - v*dt).^2 / (4*D*dt) );
    end

    function p = p3(q,dt,q0,v,D)
        p = -v/(2*D) *exp(v*q/D) .* erfc( (q + q0 + v*dt)/sqrt(4*D*dt) );
    end

    function p = p1p2p3(q,dt,q0,v,D)        
        p = p1(q,dt,q0,v,D) + p2(q,dt,q0,v,D) + p3(q,dt,q0,v,D);
    end

    function p = p2p3(q,dt,q0,v,D)        
        p = p2(q,dt,q0,v,D) + p3(q,dt,q0,v,D);
        u = 1/2 * erfc( (q0 + v*dt) / sqrt(4*D*dt) );
        p = p/u;
    end

    function p = pre(q,dt,q0,v,D)
        p = p1(q,dt,q0,v,D);
        
        if(length(q)==1)
           if (q==q0)
               p = p + 1; % supposed to be a delta function
           end
        else
           if(q(end)>=q0)
                [~,minpos] = min(abs(q-q0));
                p(minpos) = p(minpos) + 1;
           end
        end
    end

    % these should give the same values
    testVal = (q0 + v*dt) / sqrt(2*D*dt);
    pWall = normcdf( - testVal );    
    pWall2 = 0.5 * erfc( (q0 + v*dt) / sqrt(4*D*dt) );

    qRej = HSrej(dt,q0,v,D,N);

    bins = 0:0.01:2;
    binWidth = bins(2)-bins(1);
    
    countsRej = histc(qRej,bins);
    probSampledRej = countsRej/N;
    probSampledRej = probSampledRej / binWidth;
    
    figure
    bar(bins+binWidth,probSampledRej);    
    hold on

    dx = 0.0001;
    qPlot = 0:dx:2;
        
    
    p123 = p1p2p3(qPlot,dt,q0,v,D);
    
    pRej = pre(qPlot,dt,q0,v,D);
    
    %sum(f)*dx
    
    %plot(qPlot,p123,'k');
    hold on
    plot(qPlot,pRej,'r');
    
    
%     qNew = zeros(N,1);
% 
%     for j = 1:N
%         qNew(j) = newHSTest(dt,q0,v,D);
%     end

%     qNew = newHSTestN(dt,q0*ones(N,1),v,D);

    qNew = HSnew(dt,q0,v,D,N);

    countsNew = histc(qNew,bins);
    probSampledNew = countsNew/N;
    probSampledNew = probSampledNew / binWidth;
    
    figure
    bar(bins+binWidth,probSampledNew); 
    
    hold on
    
    %p23 = p2p3(qPlot,dt,q0,v,D);   
    %plot(qPlot,p23,'r');
    
    plot(qPlot,p123,'r');
    
    
end

