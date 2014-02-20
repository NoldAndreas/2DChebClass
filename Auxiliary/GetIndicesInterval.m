function I = GetIndicesInterval(x)

N = length(x);

left  = (x == x(1));
right = (x == x(end));

bound = (left | right);

nLeft  = zeros(N);
nRight = zeros(N);

nLeft(left,left) = -1;
nRight(right,right) = 1;

I = struct('left',left,'right',right,'bound',bound, ...
           'normalLeft',nLeft, ...
           'normalRight',nRight, ...
           'normal', nLeft(bound,:) + nRight(bound,:));
       
end