function C = InvMatr(A)
    [p,h1s,h2s] = size(A);
    C = zeros(p,2,2);
    
    dets     = A(:,1,1).*A(:,2,2)-A(:,1,2).*A(:,2,1);
    
    C(:,1,1) = A(:,2,2)./dets;
    C(:,1,2) = -A(:,1,2)./dets;
    C(:,2,1) = -A(:,2,1)./dets;
    C(:,2,2) = A(:,1,1)./dets;
end