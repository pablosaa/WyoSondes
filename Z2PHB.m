function PHB = Z2PHB(Z, lat, lon, alt)
	% Function to calculate the Geopotential from Geometric altitudes
	% for a given location coordinates (lat, lon, alt).
	%
	% INPUT:
	% * lat: latidude [deg]
	% * lon: longitude [deg]
	% * alt: altitude [m]
	%
	% OUTPUT:
	% * PHB: Geopotential altitude [m]
	%
	% (c) 2018, Pablo Saavedra G.
	% Geophysical Institude, University of Bergen
	% See LICENSE
	% -----------------------------------------------------------------------------

	% Defining constants:
	a = 6378137;       % [m] semi-major axes, radius on the equator
	b = 6356752.3142;  % [m] semi-minor axes, radius on the poles
	R = 6371008.8;     % [m] average radius of the sphere Earth at about 35deg latitude
	ge= 9.7803253359;  % [m/s2] effective acceleration on the equator
	gp= 9.8321849378;  % [m/s2] effective acceleration on the poles
	GM= 3.986004418e14; % [m3/s2] Geocentric gravitational constant G*M
	T = 86164098903691; % [s] rotational period (sideral day)

	omg = 2*pi/T;      % [rad/s] angular speed
	f = (a-b)/a;
	m = (b*(a*omg)^2)/GM;
	ee = (a^2 - b^2)/a^2;
	Nphi = a/sqrt(1 - ee*sind(lat).^2 );

	Z = Z(:) + alt;
	
	% Position on the Reference Ellipsoid
	Pr = (Nphi + Z)*cosd(lat);
	
	% Effective Gravitational Acceleration on the Surface
	% with g0 = effective acceleration at latitude lat
	g0 = (a*ge*cosd(lat).^2 + b*gp*sind(lat).^2)./sqrt(a*a*cosd(lat).^2 + b*b*sind(lat).^2);

	% Effective Gravitational Acceleration at Altitude
	f0 = 1;
	f1 = -2*(1 + f + m -2*f*sind(lat).^2)/a;
	f2 = 3/a^2;
	gh = @(h) g0.*(f0 + f1*h + f2*h.^2);


	ghC_vec = -gh(Z)*[cosd(lat)*cosd(lon) cosd(lat)*sind(lon) sind(lat)];
	ahC_vec = omg^2*Pr*[cosd(lon) sind(lon) 0];
	gh_vec = ghC_vec - ahC_vec;

	gZ = sqrt(gh_vec(:,1).^2 + gh_vec(:,2).^2 + gh_vec(:,3).^2);
	% Finding geopotential heights
	delZ = diff(sort([0; Z(:)]));
	
	PHB = cumsum(gZ.*delZ)/g0;
	
end
% end of function.
