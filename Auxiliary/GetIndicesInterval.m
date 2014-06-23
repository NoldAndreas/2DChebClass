function I = GetIndicesInterval(this)

x = this.Pts.x;
y = this.Pts.y;

N = length(x);

left  = (x == x(1));
right = (x == x(end));

bound = (left | right);

leftFinite  = isfinite(y(left));
rightFinite = isfinite(y(right));

finite   = ( leftFinite & left) | ( rightFinite & right);
infinite = (~leftFinite & left) | (~rightFinite & right);

nLeft  = zeros(N);
nRight = zeros(N);

nLeft(left,left) = -1;
nRight(right,right) = 1;

nFinite   = leftFinite*nLeft(finite,:) + rightFinite*nRight(finite,:);
nInfinite = (~leftFinite)*nLeft(infinite,:) + (~rightFinite)*nRight(infinite,:);

I = struct('left',left,'right',right,'bound',bound, ...
           'normalLeft',nLeft(left,:), ...
           'normalRight',nRight(right,:), ...
           'normal', nLeft(bound,:) + nRight(bound,:), ...
           'finite',finite,'infinite',infinite, ...
           'normalFinite',nFinite,'normalInfinite',nInfinite);
end