function[vertices]=gps2xyz(vertices,ellip)

%GEODETIC2ECEF Transform geodetic to geocentric (ECEF) coordinates
if ~isempty(ellip.Code)
    lon = vertices(:,1)*pi/180;
    lat = vertices(:,2)*pi/180;
    if size(vertices,2)==3
        h   = vertices(:,3);
    else
        h = lon*0;
    end
    
    % Ellipsoid has the form: [semimajor_axis, eccentricity].
    a  = ellip.SemimajorAxis;
    f  = ellip.Flattening;
    
    % Inputs phi and lambda are in radians.
    inDegrees = false;
    [rho, z] = map.geodesy.internal.geodetic2cylindrical(lat, h, a, f, inDegrees);
    [x, y]   = pol2cart(lon,rho);
    vertices = [x,y,z];
end

% This is what gps2xyz does
% N  = a/sqrt(1-e^2*sind(lat)^2);
% x  = (N+h)*cosd(lat)*cosd(lon);
% y  = (N+h)*cosd(lat)*sind(lon);
% z  = ((1-e^2)*N+h)*sind(lat);
% vertices = [x y z];