function vel = velocity_interp(xxq, yyq, ttq ,xx, yy, tt, uu, vv)
% function [uuq, vvq] = velocity_interp(xxq, yyq, t ,xx, yy, uu, vv)
    % interpolation  sur meshgrid
    
    uuq = interpn(xx,yy,tt,uu,xxq,yyq,ttq, 'linear', 0);
    vvq = interpn(xx,yy,tt,vv,xxq,yyq,ttq, 'linear', 0);
    
    vel = zeros(size(uuq,1), size(uuq,2), 2);
    vel(:,:,1) = squeeze(uuq);
    vel(:,:,2) = squeeze(vvq);
%     uuq = interp2(xx,yy,uu,xxq,yyq, 'linear', 0);
%     vvq = interp2(xx,yy,vv,xxq,yyq, 'linear', 0);
%     vel = zeros(size(uuq,1), size(uuq,2), 2);
%     vel(:,:,1) = uuq;
%     vel(:,:,2) = vvq;
end
