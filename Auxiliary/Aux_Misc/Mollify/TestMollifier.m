%function TestMollifier

    %load('/Users/NoldAndreas/Desktop/DDFT_Planar_HI-0.01-0.02-0.01-3-1-2-NaN-1-1-0-500000-1-50000-10-1-1.mat');
    y1s = -5:1:5;
    y2s = -5:1:5;

    y1 = kron(y1s,ones(size(y2s)));
    y2 = kron(ones(size(y1s)),y2s);
    tic
    z =  ApplyMollifier(x,y1,y2,@Mollifier);
    toc
    %mesh(reshape(y1',9,9),reshape(y2',9,9),reshape(z',9,9));
    
    
    mesh(reshape(y1',11,11),reshape(y2',11,11),reshape(z',11,11));
%end
