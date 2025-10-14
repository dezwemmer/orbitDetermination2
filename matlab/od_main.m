%% Orbit Determination Main
format long;

%% Function Definitions %%

% Function to calculate the DCM between ECEF and ECI
function dcm = dcmEcef2Eci(locSidTime)
    dcm = [cosd(locSidTime) -sind(locSidTime) 0;
           sind(locSidTime)  cosd(locSidTime) 0;
           0                 0                1];
endfunction


%% Variable Definitions %%

% Known Constants
Re = 6378.1366;                % Radius of Earth (km)
Ee = 0.0818191908426215;       % Eccentricity of Earth (unitless)
Omegae = 72.92115E-06;         % Angular Velocity of Earth (rad/s)
JD_prop = 2454873.2055555555; % Julian date to propagate to

% Observation Site Information
lat = 32.2227;      % degrees (geodetic)
lon = -110.0101;    % degrees
alt = 0.757;        % km

% I-Satellite Observations
JD_I = [2454872.241766892; 2454872.241940503; 2454872.242114115];
ra_I = [30.859159090717; 14.564451739639; 0.829762748762];
dec_I = [79.318796817875; 78.120651560859; 75.903618501209];
lst_I = [295.996368384251; 296.059039499724; 296.121710615405];

% C-Satellite Observations
JD_C = [2454871.514010361; 2454871.514183972; 2454871.514357583];
ra_C = [5.931355414284; 6.369337583606; 6.814572192903];
dec_C = [-26.399712354399; -23.712111094605; -20.736087662850];
lst_C = [33.286705754698; 33.349376869963; 33.412047985644];


%% Calculate Site ECEF coordinates %%

% Calculate auxiliary terms
auxC = Re / sqrt(1 - (Ee^2) * sind(lat)^2)
auxS = auxC * (1 - Ee^2)

% Calculate vertical & horizontal components
rDelta = (auxC + alt)*cosd(lat)
rK = (auxS + alt)*sind(lat)

% Calculate ECEF site position vector
rSiteEcef = [rDelta * cosd(lon);
             rDelta * sind(lon);
             rK]

% Calculate ECI site position vectors at each LST
rSiteEci = [];
for i = 1:3
  rSiteEci = [rSiteEci dcmEcef2Eci(lst_I(i)) * rSiteEcef];
endfor

rSiteEci = rSiteEci';

%% Common Computations %%

% Compute Line-Of-Site unit vectors (Luv) to the satellite at each observation time
% TODO: make this smart and modular
Luv1_I = [cosd(dec_I(1,1))*cosd(ra_I(1,1));
        cosd(dec_I(1,1))*sind(ra_I(1,1));
        sind(dec_I(1,1))];
Luv2_I = [cosd(dec_I(2,1))*cosd(ra_I(2,1));
        cosd(dec_I(2,1))*sind(ra_I(2,1));
        sind(dec_I(2,1))];
Luv3_I = [cosd(dec_I(3,1))*cosd(ra_I(3,1));
        cosd(dec_I(3,1))*sind(ra_I(3,1));
        sind(dec_I(3,1))];


%% Gauss's Technique (Angles Only)l
tau1 = JD_I(1) - JD_I(2);
tau3 = JD_I(3) - JD_I(2);

a1 = tau3 / (tau3 - tau1);
a3 = -tau1 / (tau3 - tau1);

a1u = (tau3*((tau3 - tau1)^2 - tau3^2)) / (6*(tau3 - tau1));
a3u = -(tau1*((tau3 - tau1)^2 - tau1^2)) / (6*(tau3 - tau1));

% Does this invert the same using Cramer's rule? (no)
L = [Luv1_I, Luv2_I, Luv3_I]

% Use Cramer's rule to invert
Linv = [Luv1_I(2)*Luv3_I(3)-Luv3_I(2)*Luv2_I(3)  -Luv1_I(2)*Luv3_I(3)+Luv3_I(2)*Luv1_I(3)   Luv1_I(2)*Luv2_I(3)-Luv2_I(2)*Luv1_I(3);
        -Luv2_I(1)*Luv3_I(3)+Luv3_I(1)*Luv2_I(3)  Luv1_I(1)*Luv3_I(3)-Luv3_I(1)*Luv1_I(3)  -Luv1_I(1)*Luv2_I(3)+Luv2_I(1)*Luv1_I(3);
        Luv2_I(1)*Luv3_I(2)-Luv3_I(1)*Luv2_I(2)  -Luv1_I(1)*Luv3_I(2)+Luv3_I(1)*Luv1_I(2)  Luv1_I(1)*Luv2_I(2)-Luv2_I(1)*Luv1_I(2)]

M = Linv*rSiteEci

d1 = M(2,1)*a1 - M(2,2) + M(2,3)*a3
d2 = M(2,1)*a1u + M(2,3)*a3u

C = dot(Luv2_I, rSiteEci(2,:))

% Solve for correct real root

%% Mostly Laplace -- bring this in if needed

% Compute Line-Of-Site unit vectors at any time & derivatives (L, Ld, Ldd)
% NOTE: We make the assumption that the middle time is zero (t = t2 = 0)
t = 0;
t1 = JD_I(1,1);
t2 = 0;
t3 = JD_I(3,1);
L = [((t - t2)*(t - t3)*Luv1_I)/((t1 - t2)*(t1 - t3)) ...
      + ((t - t1)*(t - t3)*Luv2_I)/((t2 - t1)*(t2 - t3))...
      + ((t - t1)*(t - t2)*Luv3_I)/((t3 - t1)*(t3 - t2))];

Ld = [((2*t - t2 - t3)*Luv1_I)/((t1 - t2)*(t1 - t3)) ...
      + ((2*t - t1 - t3)*Luv2_I)/((t2 - t1)*(t2 - t3)) ...
      + ((2*t - t1 - t2)*Luv3_I)/((t3 - t1)*(t3 - t2))];

Ldd = [(2*Luv1_I)/((t1 - t2)*(t1 - t3)) ...
      + (2*Luv2_I)/((t2 - t1)*(t2 - t3)) ...
      + (2*Luv3_I)/((t3 - t1)*(t3 - t2))];

% Compute Observer Site Accelerations
##rSited = cross([0, 0, Omegae],rSite);
##rSitedd = cross([0, 0, Omegae],rSited)


