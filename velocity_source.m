function vel = velocity_source(x,t,positions_sources, q_sources)
    %rentre un vecteur x de taille 2*n avec n le nombre de points
    
    vel = zeros(2,1);
    delta = 1;

    for i=1:length(positions_sources)
        q = q_sources{i};
        temp = positions_sources{i};
        xc = temp(1);
        yc = temp(2);

        vel(1,:) = vel(1,:) + q*(x(1,:)-xc)./(2*pi*((x(1,:)-xc).^2 + (x(2,:)-yc).^2 + delta));
        vel(2,:) = vel(2,:) + q*(x(2,:)-yc)./(2*pi*((x(1,:)-xc).^2 + (x(2,:)-yc).^2 + delta));
    end





end

