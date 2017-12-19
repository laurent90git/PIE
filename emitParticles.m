function [ xxf, yyf ] = emitParticles( xx_start,yy_start, velocityFunc, tspan )
%EMITPARTICLES computes the trajectory of particles emitted
% at all the points described in xx and yy (meshgrid-like matrices)
% and gives back xxf and yyf the positions reached at the end of tspan.
% the particles starting at xx_start(i,j) and yy_start(i,j) ends up at
% xxf(i,j) and yyf(i,j).
% Tspan, time vector in seconds
% velocityFunc gives the velocity vector for a given point.

xxf = zeros(size(xx_start));
yyf = zeros(size(xx_start));
for i=1:size(xx_start,1)
    fprintf('%.2f %%\n',i/size(xx_start,1)*100);
%     for j=1:size(xx_start,2)        
    parfor j=1:size(xx_start,2)        
        x0=[xx_start(i,j),yy_start(i,j)]';
        [t,x] = ode45(velocityFunc,tspan,x0);
%         figure, plot(x(:,1), x(:,2));
        xxf(i,j) = x(end,1);
        yyf(i,j) = x(end,2);
%         figure, plot(x(:,1),x(:,2))
    end
end


end

