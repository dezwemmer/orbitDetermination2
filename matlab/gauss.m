Re = 6378.1366;                     %km
Ee = 0.0818191908426215;
lat = 32.2227;                      %degrees
lon = -110.0101;                    %degrees
alt = 0.757;                        %km
JD_Prop = 2454873.2055555555;       %Final Julian Date to Propagate to

%Determine r_site ecef
rho = [(Re+alt)*cosd(lat)*cosd(lon) (Re+alt)*cosd(lat)*sind(lon) (Re+alt)*sind(lon)];
rho = rho'/1000;

%% SITE-TRACK
alt = 2.187 % km
lat = 39.007;
lon = -104.883;
auxC = Re / sqrt(1 - (Ee^2)*sind(lat)^2);
auxS = auxC * (1 - Ee^2);

rDelta = (auxC + alt) * cosd(39.007);
rK = (auxS + alt) * sind(39.007);

rsite = [rDelta * cosd(lon);
         rDelta * sind(lon);
         rK]
