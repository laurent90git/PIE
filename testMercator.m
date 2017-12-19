sPath = 'global-analysis-forecast-phy-001-024-hourly-t-u-v-ssh_1511272611247.nc';
ncdisp(sPath)

info = ncinfo(sPath);

data=struct();
varNames = {info.Variables.Name};
for i=1:length(varNames)
    current_field = varNames{i};
    data.(current_field) = ncread(sPath, current_field);
end

% convert to double
dataFields = fieldnames(data);
for i=1:length(dataFields)
    current_field = dataFields{i};
    data.(current_field) = double( data.(current_field) );
end

x=data.longitude;
y=data.latitude;
[xx,yy] = meshgrid(data.longitude, data.latitude);
% data.time est en heures depuis 01/01/1950

%% Plots velocity field


figure,
i=0;
while true
    i = 1+mod(i,length(data.time));
    % for i=1:length(data.time)
    %    quiver(xx',yy',data.uo(:,:,:,i), data.vo(:,:,:,i))
    surf(xx',yy',(data.uo(:,:,:,i).^2 + data.vo(:,:,:,i).^2).^0.5)
    shading flat;
    view([0 0 1])
    drawnow;
    %    pause(0.2);
end

%% Emit particles
x_start = data.longitude(180:5:200);
y_start = data.latitude(180:5:200);
[xx_start, yy_start] = meshgrid(x_start, y_start);
% start positions are latitude and longitude, hence we need to convert the
% speed to deg/s of latitude and longitude
tspan = data.time - min(data.time);

%% convert speed to deg/h (coefficient depends on latitude)
r_earth = 6378137; %meters
factor_convert_angular_speed_x = 1./(r_earth*cos(pi/180*yy' ));
factor_convert_angular_speed_y = 1./(r_earth);

uu = 3600 * data.uo .* repmat(factor_convert_angular_speed_x, 1,1,size(data.uo,3), size(data.uo,4));
vv = 3600 * data.vo .* repmat(factor_convert_angular_speed_y,size(data.vo,1),size(data.vo,2),size(data.vo,3), size(data.vo,4));

% surf(xx',yy',(uu(:,:,:,i).^2 + vv(:,:,:,i).^2).^0.5); shading flat;
[xxx,yyy,ttt] = ndgrid(x, y, tspan);

% figure, surf(xx',yy',factor_convert_angular_speed_x); shading flat;
%% Time integration
velocityFunc = @(t,x) squeeze(velocity_interp(x(1,:), x(2,:), t ,xxx, yyy, ttt, squeeze(uu), squeeze(vv)));
[ xxf, yyf ] = emitParticles( xx_start,yy_start, velocityFunc, tspan );

velocity_interp(xx_start(1,:), yy_start(1,:), 0 ,xxx, yyy, ttt, squeeze(uu), squeeze(vv))
% velocityFunc(0, [xx_start(1,1), yy_start(1,1)]);

distances_parcourues = r_earth*(xx_start-xxf)/1e3; % d en km
%% Plots
figure,
for ii=1:size(xx_start,1)
    for jj=1:size(xx_start,2)
        hold on;
        scatter(xx_start(ii,jj), yy_start(ii,jj), 'o');
        scatter(xxf(ii,jj), yyf(ii,jj), '+');
%         plot([xx_start(ii,jj), yy_start(ii,jj)], [xxf(ii,jj), yyf(ii,jj)]);
    end
end