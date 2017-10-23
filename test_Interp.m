x = 1:10;
y = x;

[xx,yy] = meshgrid(x,y);
zz = (xx-5).^2 + (yy-5).^2;
uu = zz;
vv= -zz;
t = 0;

% xq = 2:0.25:6;
% yq = 2:0.1:7;
% [xxq,yyq] = meshgrid(xq,yq)
% 
% zzq = interp2(xx,yy,zz,xxq,yyq)
% 
% 
% figure,
% subplot(2,1,1)
% surf(xx,yy,zz)
% hold on
% surf(xxq,yyq,zzq)


xxq = 2.63;
yyq = 8.25;
vel = squeeze(velocity_interp(xxq, yyq, t ,xx, yy, uu, vv))