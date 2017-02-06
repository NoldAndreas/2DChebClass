function plotParticlePaths2D(x)

nParticles = size(x,1)/2;

hf = figure;
ha = axes;

for iParticle = 1:nParticles
    xi = x((iParticle-1)*2 +1,:);
    yi = x(iParticle*2,:);
    
    plot(ha,xi(1),yi(1),'o');
    hold on
    plot(ha,xi(end),yi(end),'x');
    plot(ha,xi,yi,'-')
    
end