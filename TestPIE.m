% Je sais Laurent, c'est pas optimisé mais c'est un premier code. Soit
% indulgent


%% Petit test pour voir les champs de vitesse

clear all;  
% close all;

x = [-1:0.02:1];
y = [-1:0.02:1];
[xx,yy] = meshgrid(x,y);

% positions_vortices = {[-0.5,-0.5],[-0.5,0.5],[0.5,0.5],[0.5,-0.5]};
% gamma_vortices  =    {1,-1,1,-1};
positions_vortices = {[0,0]};
gamma_vortices  =    {-1};


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


tspan = [0 10];
x0 = [0 0.2];

velocityFunc = @(t,x) squeeze(velocity_interp(x(1), x(2), t ,xx, yy, squeeze(velocityField(:,:,1)), squeeze(velocityField(:,:,2))));
% velocityFunc = @(t,x) velocity_interp(x,t,xx, yy,velocityField);

[t,x] = ode45(velocityFunc,tspan,x0);

figure;
plot(x(:,1),x(:,2));

%% calcul des surfaces cohÃ©rentes par le mÃ©thode "LAVD"

%Implémenté uniquement pour les vortex et pour un rotationnel moyen nul
%(deux vortexs d'intensité opposées par exemple). J'ai pas l'impression que
%ce soit une méthode très utile pour nous.

% tspan = [0 10];
% [xx,yy] = meshgrid([-1:0.05:1,-1:0.05:1]);
% dx = 0.01;
% dy = 0.01;
% LAVD = zeros(size(xx,1),size(xx,2));
% 
% for i=1:length(xx)
%     i/length(xx)*100
%     for j=1:length(yy)
%         x0=[xx(i,j),yy(i,j)];
%         [t,x] = ode45(velocityFunc,tspan,x0);
%         x = x';
%         
%         dvdx = (velocity_vortex(x+dx*[1;0]*ones(1,size(x,2)),t,positions_vortices, gamma_vortices) - velocity_vortex(x,t,positions_vortices, gamma_vortices))/dx;
%         dudy = (velocity_vortex(x+dy*[0;1]*ones(1,size(x,2)),t,positions_vortices, gamma_vortices) - velocity_vortex(x,t,positions_vortices, gamma_vortices))/dy;
%         omega = dvdx - dudy; %rotationel pour un point initial pour tout les temps
%         
%         for k=1:length(t)-1
%             LAVD(i,j) = LAVD(i,j) + (norm(omega(:,k)) -1)*(t(k+1)-t(k));
%         end
%         
%     end 
% end
% 
% figure;
% contourf(xx,yy,LAVD);
% surf(xx,yy,LAVD);
% colorbar;

%% Cauchy-Green:
%Méthode qui donne des résultats pas trop idiots

tspan = [0 5];
xi1 = zeros(size(xx,1),size(xx,2),2);
xi2 = zeros(size(xx,1),size(xx,2),2);
ev1 =  zeros(size(xx,1),size(xx,2));
ev2 =  zeros(size(xx,1),size(xx,2));
scalar =  zeros(size(xx,1),size(xx,2));

% dx = 0.01;
% dy = 0.01;


xxf = zeros(size(xx));
yyf = zeros(size(xx));
for i=1:size(xx,1)
    i/size(xx,1)*100
    parfor j=1:size(xx,2)        
        x0=[xx(i,j),yy(i,j)]';
        [t,x] = ode45(velocityFunc,tspan,x0);
        xxf(i,j) = x(end,1);
        yyf(i,j) = x(end,2);
    end
end

%%
figure,
surf(xx,yy,yyf)
%% 
[dxF_x, dxF_y] = gradient(xxf, xx(1,:), yy(:,1));
[dyF_x, dyF_y] = gradient(yyf, xx(1,:), yy(:,1));


for i=1:size(xx,1)
    i/size(xx,1)*100
    for j=1:size(xx,2)   
        deltaF = [dxF_x(i,j) , dxF_y(i,j)  ;
                  dyF_x(i,j) , dyF_y(i,j) ];
        
        C = deltaF'*deltaF;
        [V,D] = eig(C);
        if D(1,1)>D(2,2)
            xi1(i,j,:) = V(:,2);
            xi2(i,j,:) = V(:,1);
            ev1(i,j) = D(2,2);
            ev2(i,j) = D(1,1);
            
        else
            xi1(i,j,:) = V(:,1);
            xi2(i,j,:) = V(:,2);
            ev1(i,j) = D(1,1);
            ev2(i,j) = D(2,2);
        end
    end
end

%%
%         deltaF = [dx1dx0 dx1dy0;
%             dy1dx0 dy1dy0];
%         
%         C = deltaF'*deltaF;
%         [V,D] = eig(C);
%         if D(1,1)>D(2,2)
%            xi1(i,j,:) = V(:,2);
%            xi2(i,j,:) = V(:,1);
%            ev1(i,j) = D(2,2);
%            ev2(i,j) = D(1,1);
%            
%         else
%            xi1(i,j,:) = V(:,1);
%            xi2(i,j,:) = V(:,2);
%            ev1(i,j) = D(1,1);
%            ev2(i,j) = D(2,2);
%         end  


% for i=1:length(xx)
%     i/length(xx)*100
%     for j=1:length(yy)
%         x0=[xx(i,j),yy(i,j)];
%         %[t,x] = ode45(odefun,tspan,x0);
%         [tdx,xdx] = ode45(velocityFunc,tspan,x0+dx*[1 0]);
%         [tmdx,xmdx] = ode45(velocityFunc,tspan,x0-dx*[1 0]);
%         [tdy,xdy] = ode45(velocityFunc,tspan,x0+dy*[0 1]);
%         [tmdy,xmdy] = ode45(velocityFunc,tspan,x0-dy*[0 1]);
%                 
%         %deformation gradient tensor
%         dx1dx0 = (xdx(size(xdx,1),1)- xmdx(size(xmdx,1),1))/(2*dx);
%         dx1dy0 = (xdy(size(xdy,1),1)- xmdy(size(xmdy,1),1))/(2*dy);
%         dy1dx0 = (xdx(size(xdx,1),2)- xmdx(size(xmdx,1),2))/(2*dx);
%         dy1dy0 = (xdy(size(xdy,1),2)- xmdy(size(xmdy,1),2))/(2*dy);
%         
%         deltaF = [dx1dx0 dx1dy0;
%             dy1dx0 dy1dy0];
%         
%         C = deltaF'*deltaF;
%         [V,D] = eig(C);
% 
%         if D(1,1)>D(2,2)
%            xi1(i,j,:) = V(:,2);
%            xi2(i,j,:) = V(:,1);
%            ev1(i,j) = D(2,2);
%            ev2(i,j) = D(1,1);
%            
%         else
%            xi1(i,j,:) = V(:,1);
%            xi2(i,j,:) = V(:,2);
%            ev1(i,j) = D(1,1);
%            ev2(i,j) = D(2,2);
%         end  
%         scalar(i,j) = abs(velocity_vortex(x(size(x,1),:)',t,positions_vortices, gamma_vortices)'*(deltaF*[xi1(i,j,1); xi1(i,j,1)]));
%     end 
% end

etaPlus = zeros(size(xx,1),size(xx,2),2);
etaPlus(:,:,1) = ((ev2 -1)./(ev2-ev1)).^(0.5).*xi1(:,:,1) + ((1 - ev1)./(ev2-ev1)).^(0.5).*xi2(:,:,1);
etaPlus(:,:,2) = ((ev2 -1)./(ev2-ev1)).^(0.5).*xi1(:,:,2) + ((1 - ev1)./(ev2-ev1)).^(0.5).*xi2(:,:,2);
etaPlus = real(etaPlus);

etaMoins = zeros(size(xx,1),size(xx,2),2);
etaMoins(:,:,1) = ((ev2 -1)./(ev2-ev1)).^(0.5).*xi1(:,:,1) - ((1 - ev1)./(ev2-ev1)).^(0.5).*xi2(:,:,1);
etaMoins(:,:,2) = ((ev2 -1)./(ev2-ev1)).^(0.5).*xi1(:,:,2) - ((1 - ev1)./(ev2-ev1)).^(0.5).*xi2(:,:,2);
etaMoins = real(etaMoins);
%% DIfférents plots des résultats
% 
% figure;
% quiver(xx0,yy0,xi1(:,:,1),xi1(:,:,2),'b');
% legend('\xi_1');
% 
figure;
quiver(xx,yy,xi2(:,:,1),xi2(:,:,2),'r');
legend('\xi_2')


% figure;
% quiver(xx,yy,etaPlus(:,:,1),etaPlus(:,:,2),'r')
% legend('\eta_{+}');
% 
% figure;
% quiver(xx,yy,etaMoins(:,:,1),etaMoins(:,:,2),'b')
% legend('\eta_{-}')

figure;
%contourf(xx,yy,ev1);
colorbar;
surf(xx,yy,ev1);
axis([-0.5 0.5 -0.5 0.5 0 1]);

% 
% xlabel('x-axis');
% ylabel('y-axis');

%% Prendre une structure cohérente
10

sspan = [0 0.5];
x0 = [0.5; 0.1];

%odefun = @(t,x) velocity_vortex(x,t,positions_vortices, gamma_vortices);
% coherentStructVel = @(t,x) velocity_interp(x,t,xx,yy,etaPlus);
coherentStructVel = @(t,x) squeeze(velocity_interp(x(1), x(2), t ,xx, yy, squeeze(xi2(:,:,1)), squeeze(xi2(:,:,2))));

[sStruct,xStruct] = ode45(coherentStructVel,sspan,x0);

figure;
plot(xStruct(:,1),xStruct(:,2));
axis([-1 1 -1 1])

%% Convection de la structure

figure;


dt = 0.2;
tmax =10;
N = floor(tmax/dt);
xStructUp = xStruct';


for i=1:N
    fprintf('\n=== i/N*100 %')
    parfor j=1:size(xStructUp,2)
        [intermt,intermx] = ode45(velocityFunc,[0 dt],xStructUp(:,j));
        xStructUp(:,j) = intermx(size(intermx,1),:)';
    end
    
    plot(xStructUp(1,:),xStructUp(2,:));
    dX = xStructUp(:,2:end) - xStructUp(:,1:end-1);
    longueur = sum(sqrt(sum(dX.^2,1)));
    fprintf('\n\t l = %d',longueur);
    axis([-1.2 1.2 -1.2 1.2])
    drawnow;
    dt*N
%     pause(0.01)
    
end





