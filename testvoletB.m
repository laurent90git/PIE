%% Petit test pour voir les champs de vitesse

clear all;  
% close all;
dx = 0.02;
dy=dx;
x = [-1:dx:1];
y = [-1:dy:1];
[xx,yy] = meshgrid(x,y);

% positions_vortices = {[-0.5,-0.5],[-0.5,0.5],[0.5,0.5],[0.5,-0.5]};
% gamma_vortices  =    {1,-1,1,-1};
positions_vortices = {[-0.5,0], [0.5,0]};
gamma_vortices  =    {-1, 1};


u = zeros(size(xx,1),size(xx,2));
v = zeros(size(xx,1),size(xx,2));

delta = 0.1;

for i=1:length(positions_vortices)
    gamma = gamma_vortices{i};
    temp = positions_vortices{i};
    xc = temp(1);
    yc = temp(2);
    
    %vortices
     u = u - gamma*(yy-yc)./(2*pi*((xx-xc).^2 + (yy-yc).^2 + delta));
     v = v + gamma*(xx-xc)./(2*pi*((xx-xc).^2 + (yy-yc).^2 + delta));
%     
% %     %source
%      u = u + gamma*(xx-xc)./(2*pi*((xx-xc).^2 + (yy-yc).^2 + delta));
%      v = v + gamma*(yy-yc)./(2*pi*((xx-xc).^2 + (yy-yc).^2 + delta));
end


% u = yy.*exp(-(xx.^2 + yy.^2));
% v = -xx.*exp(-(xx.^2 + yy.^2));
% 
figure;
quiver(xx,yy,u,v);

velocityField = zeros(size(xx,1),size(xx,2),2);
velocityField(:,:,1) = u;
velocityField(:,:,2) = v;

%% test ode45


tspan = [0 20];
x0 = [0.2 0.2];

% velocityFunc = @(t,x) squeeze(velocity_interp(x(1), x(2), t ,xx, yy, squeeze(velocityField(:,:,1)), squeeze(velocityField(:,:,2))));
velocityFunc = @(t,x) velocity_vortex(x,t,positions_vortices, gamma_vortices);

[t,x] = ode45(velocityFunc,tspan,x0);

figure;
plot(x(:,1),x(:,2));

%% Throw many particles
xxf = xx; %zeros(size(xx));
yyf = yy; %zeros(size(xx));
deltaT_plot = 5;
current_tspan = [tspan(1) deltaT_plot];

figure
while current_tspan(end)<tspan(end)
    current_tspan = current_tspan + deltaT_plot;
    for i=1:size(xx,1)
        i/size(xx,1)*100
        parfor j=1:size(xx,2)
            x0=[xxf(i,j),yyf(i,j)]';
            [t,x] = ode45(velocityFunc,current_tspan,x0);
            xxf(i,j) = x(end,1);
            yyf(i,j) = x(end,2);
        end
    end
    figure(30)
    density = zeros(size(xx));
    nTotalParticles = numel(xx);
    for i=1:size(xx,1)
        for j=1:size(xx,2)
            tempDistance = (xxf-xx(i,j)).^2 + (yyf-yy(i,j)).^2 - (dx^2 + dy^2);
            temp2 = tempDistance<0;
            density(i,j) = sum(sum(temp2))/nTotalParticles;
        end
    end
    % surf(xx,yy,density);
    contourf(xx,yy,density);
    axis square
    caxis([0 1e-3])
    colorbar
    drawnow
end

%% Plot final density
figure;
plot(xxf,yyf,'+', 'Color', 'b')

%% Density computation
density = zeros(size(xx));
nTotalParticles = numel(xx);
for i=1:size(xx,1)
    for j=1:size(xx,2)
        tempDistance = (xxf-xx(i,j)).^2 + (yyf-yy(i,j)).^2 - (dx^2 + dy^2);
        temp2 = tempDistance<0;
        density(i,j) = sum(sum(temp2))/nTotalParticles;
    end
end

figure;
% surf(xx,yy,density);
contourf(xx,yy,density);
axis square
caxis([0 1e-3])
colorbar