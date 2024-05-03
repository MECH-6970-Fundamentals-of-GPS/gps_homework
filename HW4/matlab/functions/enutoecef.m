function [C] = enutoecef(lat,lon)
%ENUTOECEF Summary of this function goes here
%   Detailed explanation goes here
    C = [-sind(lon), -cosd(lon)*sind(lat), cosd(lon)*cosd(lat);
          cosd(lon), -sind(lon)*sind(lat), sind(lon)*cosd(lat);
                  0,            cosd(lat),           sind(lat)];
end

