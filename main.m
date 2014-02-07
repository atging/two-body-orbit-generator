% Andrew Ging, ASTE-580, 11/30/2011
% Computer Project
%% Description %%
% Takes position and velocity vectors and central body information as 
% inputs. Determines type of orbit, Keplerian elements, and whether the 
% spacecraft is on an impactor trajectory. Computes the position and 
% velocity vectors at a user input specified time.
clear;
format compact;
xero = 5E-4;    % tolerance for values close to zero
oType = 0;      % Orbit type marker, initialization
it = 60;      % seconds, time increment for part 4, method 2

AU = 149597871; % km, astronomical unit

%% User Inputs %%
% mu = input('Gm of Body, km^3/s^2: ');
% r0 = input('Initial position vector, km (eg. [1,2,3]): ');
% v0 = input('Initial velocity vector, km/s (eg. [4,5,6]): ');
% Rb = input('Central body mean surface radius, km: ');
% dt = input('Delta time (for t > t0), seconds: ');

% Case 0:
% MSL's Centaur vehicle ephemeris output from Horizons on 2012-Feb-29
mu = 132712440017.987;      % km^3/s^2 (Sun)
% r0 = [-1.618950631518551E+08  8.916064793611209E+07 -5.351713299778852E+06];
% v0 = [-1.795683661278650E+01 -2.028244240873942E+01 -1.945344663229576E-01];
Rb = 696000;        % km (Sun)
dt = 219-60;        % days, delta-t (time since t = t0)
dt = dt * 3600*24;  % seconds
radiansCircle = 0:(2*pi/90):(2*pi);
circleEarthX = AU * cos(radiansCircle);
circleEarthY = AU * sin(radiansCircle);
circleEarthZ = zeros(size(circleEarthY));
Heliocentric = 1;   % tells it to plot Earth's rough orbit for reference

% Mars' Orbit (parameters as of 2011-Nov-27 00:00:00)
r0 = [-1.232877989784380E+08  2.109814671850678E+08  7.447932575332564E+06];
v0 = [-2.000124286766004E+01 -1.016348705053261E+01  2.781521225303251E-01];

% Case 1:
% mu = 398600.433;    % km^3/s^2
% r0 = [-14192.498, -16471.197, 1611.2886];      % km
% v0 = [-4.0072937, -1.2757932, 1.9314620];      % km/s
% Rb = 6378.14;       % km, central body mean surface radius
% dt = 8.0;           % hours, delta-t (time since t = t0)
% dt = dt * 3600;     % seconds

% Case 2:
% mu = 132712440017.987;      % km^3/s^2
% r0 = [148204590.0357, 250341849.5862, 72221948.8400];   % km
% v0 = [-20.5065125006, 7.8793469985, 20.0718337416];     % km/s
% Rb = 696000;        % km
% dt = 10;     % days
% dt = dt * 24 * 3600;    % seconds

% Case 3:
% mu = 37940626.1;    % km^3/s^2
% r0 = [-321601.0957, -584995.9962, -78062.5449];     % km
% v0 = [8.57101142, 7.92783797, 1.90640217];          % km/s
% Rb = 60268;         % km
% dt = (24-14)*3600 + 47*60 + 39.3;       % seconds

% Case 4:
% mu = 8978.1382;         % km^3/s^2
% r0 = [8193.2875, -21696.2925, 7298.8168];       % km
% v0 = [-2.29275936, 4.94003573, -1.67537281];    % km/s
% Rb = 2575;              % km
% dt = 3600 + 4*60 + 1.18;    % seconds

% Case 5 (ECI):
% mu = 398600.433;
% r0 = [5492.00034, 3984.00140, 2.95581];
% v0 = [-3.931046491, 5.498676921, 3.665980697];
% Rb = 6378.14;
% dt = 5.0;       % hours
% dt = dt * 3600; % sec

%% (1) Determination of orbit type %%

h0 = cross(r0, v0);    % angular momentum
h0_mag = norm(h0);     % magnitude of h0
r0_mag = norm(r0);     % magnitude of r0
v0_mag = norm(v0);     % magnitude of v0
beta0 = acosd(h0_mag/(r0_mag*v0_mag));   % Flight path angle, degrees

if ( beta0 > 90 || beta0 < -90 )
    error(['Possible problem: beta0 = ', num2str(beta0)]);
end

X0 = r0_mag * v0_mag^2 / mu;
e = sqrt( ((X0-1)^2)*(cosd(beta0))^2 + (sind(beta0))^2 );   % eccentricity

if ( e <= xero && e >= 0 )  % e is approximately Zero
    disp('(1) The orbit type is Circular');
elseif ( e > (1+xero) )     % e is greater than 1 (with tolerance)
    disp('(1) The orbit is Hyperbolic');
    oType = 3;
elseif ( e >= 1 )     % e is approximately 1
    disp('(1) The orbit is Parabolic');
    oType = 2;
elseif ( e > xero )         % e is between Zero and 1 (with tolerance)
    disp('(1) The orbit is Elliptical');
    oType = 1;
else    % error catching
    disp(['e is less than zero (e = ', num2str(e), ')']);
    error('This is a problem.');
end

%% (2) Keplerian Elements Computation (t = t0) %%
% (a, b and c) For all orbit types:

p = h0_mag^2 / mu;      % km, parameter length
rp = p / (1 + e);       % km, radius at periapsis
theta = acosd( (X0*cosd(beta0)^2 - 1) / e );    % deg, true anomaly
if ( beta0 < 0 )
    theta = -abs(theta); % beta0 < 0 indicates that -180 < theta < 0
else
    theta = abs(theta); % beta0 > 0 indicates that 0 < theta < 180
end
h0_hat = h0 / h0_mag;       % unit vector of h0
K_hat = [0, 0, 1];          % unit vector of Z-direction 
h3 = dot(h0_hat, K_hat);    % projection of h0 onto z-axis
i = abs(acosd(h3));         % deg, i is between 0 and 180 deg

% N = line of nodes vector in the X, Y plane and orbit plane
N = cross(K_hat, h0_hat);   % vector of the ascending node
N_hat = N ./ norm(N);       % unit vector of N
I_hat = [1, 0, 0];          % unit vector of X-direction
Omega = acosd(dot(N_hat, I_hat));   % deg, R.A. of ascending node
if ( N_hat(2) < 0 )     % Omega in 3rd or 4th quadrants
    Omega = 360 - Omega;
end
e_vec = mu^(-1) * ((v0_mag^2-mu/r0_mag).*r0 - dot(r0,v0).*v0);  % e vector
lomega = acosd( dot(e_vec, N_hat) / e );
if ( e_vec(3) < 0 ) 
    lomega = 360 - lomega;
end

disp('(2) Keplerian Elements at initial time t0:');
% (c) Hyperbolic orbits:

if ( oType == 3 )     % If hyperbolic
    a = -p / (e^2 - 1);          % km, semi-major axis
    % hyperbolic approach and departure velocities
    v_inf = sqrt( -mu / a );    % km/s, V infinite
    table = [e, p, a, v_inf, rp, i, Omega, lomega, theta]; % elements table
    % names of first 3 elements:
    fprintf('  %s\t\t%s\t\t\t%s\n','e (unitless)','p (km)','a (km)');
    disp(['  ',num2str(table(1:3))]);  % display elements table in cmd window
    % names of last 6 elements:
    fprintf('  %s\t%s\t\t\t%s\t\t%s\t\t%s\t\t%s\n',...
        'V_infinite (km/s)','rp (km)','i (deg)','Omega (deg)',...
        'omega (deg)','theta (deg)');
    disp(['  ',num2str(table(4:9))]);  % display last 6 elements
    
% (a) Elliptical orbits:

elseif ( oType <= 1 )     % If elliptical
    a = p / (1 - e^2);              % km, semi-major axis
    tau = 2*pi*sqrt( a^3 / mu );    % seconds, orbital period
    tau_min = tau / 60;             % minutes
    ra = 2 * a - rp;                % km, radius at apoapsis
    % Eccentric anomaly:
    E0 = 2*atan2( (sqrt(1-e)*tand(theta/2)), sqrt(1+e) );   % rad, E at t0
    tp = tau_min/(2*pi)*(E0 - e*sin(E0));   % minutes, time since periapsis
    M0 = 180/pi * (2*pi / tau_min) * tp;    % deg, Mean anomaly at t0
    table = [e, p, a, tau_min, rp, ra, i, Omega, lomega, theta, ...
        180/pi*E0, M0, tp];   % elements table for display output
    % names of first 6 elements:
    fprintf('  %s\t\t%s\t\t\t%s\t\t%s\t%s\t\t\t%s\n',...
        'e (unitless)','p (km)','a (km)','tau (minutes)','rp (km)',...
        'ra (km)');
    disp(['  ',num2str(table(1:6))]);  % display first 6 elements
    % names of last 7 elements:
    fprintf('  %s\t\t%s\t\t%s\t%s\t\t%s\t\t\t%s\t\t%s\n',...
        'i (deg)','Omega (deg)','omega (deg)','theta (deg)',...
        'E (deg)','M (deg)','tp (minutes)');
    disp(['  ',num2str(table(7:13))]); % display last 7 elements
    
% (b) Parabolic orbits:
    
elseif ( oType == 2 )
    table = [p, rp, i, Omega, lomega, theta];   % elements table
    % names of elements:
    fprintf('  %s\t\t\t%s\t\t\t%s\t%s\t%s\t%s\n','p (km)',...
        'rp (km)','i (deg)','Omega (deg)','omega (deg)','theta (deg)');
    % display elements in cmd window:
    fprintf('  %d\t%d\t%3.4f\t%3.4f\t\t%3.4f\t\t%f\n',p,rp,i,Omega,lomega,...
        theta);    
end

%% (3) Spacecraft is on an impactor? %%

if ( rp <= Rb )
    disp('(3) Spacecraft WILL impact the central body ');
    disp(['  since r_p = ', num2str(rp), ...
        ' km, which is less than R_b = ', num2str(Rb), ' km']);
else
    disp('(3) Spacecraft will NOT impact the central body');
end

%% (4) Computation of position & velocity vectors at specified t > t0 %%

%% Two methods: a) 1st method; Determine (mostly) analytically:
% t-tp => E (by iteration) => theta => r_orbital => r_inertial (elliptical)
%                                      v_orbital => v_inertial also
%dt = tau;

disp('(4) Position vector *r = *r(t) & velocity vector *v = *v(t)');
disp('    at specified time t, where t > t0');

if ( oType == 3 )   % Hyperbolic
    %(not reported because it doesn't produce correct results in all cases)
    % Hyperbolic "eccentric" anomaly:
    sinhF = sqrt(e^2 - 1) * sind(theta) / (e*cosd(theta) + 1);  % sinh of F
    F0 = 2*atanh( sqrt((e-1)/(e+1))*tand(theta/2) );    % alternate F
    while ( sinhF >= 2*pi )
        sinhF = sinhF - 2*pi;   % bring down to 0 < sinhF < 2pi rad
    end
    N_h0 = e*sinhF - asinh(sinhF);  % rad, "mean" anomaly (hyperbolic)
    n = sqrt(mu/abs(a^3));              % rad/s, mean motion
    F = funcF(dt, n, N_h0, e, xero);    % rad, *my* function used
    theta_t = 2*atan2( (sqrt(e+1)*tanh(F/2)), sqrt(e-1) ); % rad, new theta

elseif ( oType <= 1 )   % Elliptical
    % rad, *my* function returns new E
    E = funcE(dt, tau, M0, e, xero, oType); 
    theta_t = 2*atan2( (sqrt(1+e)*tan(E/2)), sqrt(1-e) );  % rad, new theta
end
if ( oType == 3 || oType <= 1 ) % if not parabolic
    r_mag_n = p / (1 + e*cos(theta_t));     % km, new radius magnitude
    r_Xi = r_mag_n*cos(theta_t);            % km, r Xi component
    r_Eta = r_mag_n*sin(theta_t);           % km, Eta component
    v_Xi = sqrt(mu/p)*(-sin(theta_t));      % km/s, v Xi component
    v_Eta = sqrt(mu/p)*(e+cos(theta_t));    % km/s, v Eta component

    % alpha*T (alpha Asterisk transpose) (phi=omega, psi=Omega, theta=i)
    alphaS = ...
        [ (cosd(lomega)*cosd(Omega)-sind(lomega)*cosd(i)*sind(Omega)), ...
        (-sind(lomega)*cosd(Omega)-cosd(lomega)*cosd(i)*sind(Omega)), ...
        (sind(Omega)*sind(i));
        (cosd(lomega)*sind(Omega)+sind(lomega)*cosd(i)*cosd(Omega)), ...
        (-sind(lomega)*sind(Omega)+cosd(lomega)*cosd(i)*cosd(Omega)), ...
        (-sind(i)*cosd(Omega));
        (sind(lomega)*sind(i)), (cosd(lomega)*sind(i)), cosd(i) ];
    % position & velocity vectors after t = dt + t0:
    rJ2000 = alphaS * [r_Xi; r_Eta; 0];     % km, r vector inertial
    vJ2000 = alphaS * [v_Xi; v_Eta; 0];     % km/s, v vector inertial
    if ( oType <= 1 )
        % display the *r(t) and *v(t) vectors if elliptical
        disp(['  At t = ', num2str(dt), ' seconds > t0;']);
        disp('    Determined via method 1:');
        disp(['    (t-tp => E (by iteration) => theta => r_orbital => ',...
            'r_inertial)']);
        disp(['    *r(t) = [',num2str(transpose(rJ2000)),'] km']);
        disp('    and v_inertial:');
        disp(['    *v(t) = [',num2str(transpose(vJ2000)),'] km/s']);
    end
end

%% Verification:
% b) 2nd method; Small time increments and corresponding variations in 
% the position and velocity vectors w/time.
% 2nd verification: do methods 1 & 2 for t = tau + t0, which should return
% same position and velocity vectors as at t = t0 for elliptical orbits

% 2nd method:
% Initial loop variables
t = 0;          % seconds, initializing time variable
r = r0;         % km, position vector
v = v0;         % km/s, velocity vector
if ( (dt/3600) > 24 )   % if dt > one day
    rInt = 3600;    % seconds, record interval for plot
elseif ( (dt/3600) > 4 )
    rInt = 500;    
else
    rInt = 10;
end
numInt = floor(dt/it/rInt);     % Number of intervals for plotting
r_mat = zeros(numInt, 3);       % Initializing position matrix for plot
if ( oType <= 1 )
    if ( dt < tau ) % for plotting remainder of orbit in diff color
        r_mattau = zeros(floor((tau-dt)/it/rInt), 3);    % initialize
    end
end
%v_mat = r_mat;                  
count = 0;

while ( t < dt )
    r_mag = norm(r);            % km, magnitude of position
    v_mag = norm(v);            % km/s, magnitude of velocity
    dv_mag = mu / r_mag^2 * it; % km/s, incremental delta-v magnitude
    r_hat = -r ./ r_mag;        % unit vector of dv direction
%    v_hat = v ./ v_mag;         % unit vector of velocity
    dv = dv_mag .* r_hat;       % km/s, delta-v vector (x,y,z coords)
    v = v + dv;                 % km/s, new velocity after time increment
    dr = v * it;                % km, vector of delta-position lengths
    r = r + dr;                 % km, new position after time increment it
    if ( mod(count, rInt) == 0 && t <= dt )    % for plot of t0 < t < dt
        % record position into matrix for plot
        r_mat((count/rInt+1),1:3) = r;
        %v_mat((count/rInt+1),1:3) = v;
    end
    if ( (dt-t) <= 2*it && (dt-t) >= 0 )
        rfin = r;
        vfin = v;
    end    
    t = t + it;                 % seconds, new time (incremented)
    count = count + 1;
    
    if ( oType <= 1 )   % only for elliptical orbits is this possible
        % keep going til t=(t0+tau) for rest of plot
        if ( dt < tau && t >= dt )
            countf = count;
            while ( t < tau )   % same code as above
                r_mag = norm(r);            
                v_mag = norm(v);            
                dv_mag = mu / r_mag^2 * it; 
                r_hat = -r ./ r_mag;        
                v_hat = v ./ v_mag;         
                dv = dv_mag .* r_hat;       
                v = v + dv;                 
                dr = v * it;               
                r = r + dr;
                % for plot of t0 < t < (t0 + tau)
                if ( mod((count-countf), rInt) == 0 )    
                    % record position into matrix for plot
                    r_mattau(((count-countf)/rInt+1),1:3) = r;
                    %v_mattau(((count-countf)/rInt+1),1:3) = v;
                end
                t = t + it; 
                count = count + 1;
            end
        end
    end
    
end
disp(' ');  % display *r(t) and *v(t) determined iteratively:
disp('    Determined via method 2 (iteratively):');
disp(['    *r(t) = [',num2str(r),'] km']);
disp(['    *v(t) = [',num2str(v),'] km/s']);

%% (v) Something extra related to project: Plot orbits %%
figure(1);
% Plot orbital positions throughout t0 < t < dt
plot3(r_mat(1:numInt,1), r_mat(1:numInt,2), r_mat(1:numInt,3),'r');
hold on;
if ( oType <= 1  )   % elliptical specific stuff
    if ( dt < tau )
        Lr = length(r_mattau);
        plot3(r_mattau(1:Lr,1), r_mattau(1:Lr,2), r_mattau(1:Lr,3),'--k');
        legend('t_0 < t < dt', 'dt < t < tau');
    else
        legend('t_0 < t < dt');
    end
else
    legend('t_0 < t < dt');
end
[Xs,Ys,Zs] = sphere;    % generate a sphere for the central body
mesh(Xs.*Rb,Ys.*Rb,Zs.*Rb);     % plot the central body
grid on;
xlabel('X axis (km)');
ylabel('Y axis (km)');
zlabel('Z axis (km)');
plot3([-Rb*2 Rb*2], [0 0], [0 0],'k');  % x-axis display
plot3([0 0], [-Rb*2 Rb*2], [0 0],'k');  % y-axis display
plot3([0 0], [0 0], [-Rb*2 Rb*2],'k');  % z-axis display

if ( Heliocentric )
    plot3(circleEarthX, circleEarthY, circleEarthZ);    % plot Earth's orbit
end

%Final values (checking of hyperbolic anomaly calculations; not reported):
if ( oType == 3 )
    hf = cross(r, v);    % angular momentum
    h_magf = norm(hf);    % magnitude of h0
    r_magf = norm(r);     % magnitude of r0
    v_magf = norm(v);     % magnitude of v0
    betaf = acosd(h_magf/(r_magf*v_magf));   % Flight path angle, degrees
    Xf = r_magf * v_magf^2 / mu;
    thetaf = acosd( (Xf*cosd(betaf)^2 - 1) / e );    % deg, true anomaly
    sinhFf = sqrt(e^2 - 1) * sind(thetaf) / (e*cosd(thetaf) + 1);
    Ff = 2*atanh( sqrt((e-1)/(e+1))*tand(thetaf/2) );    % rad
    N_hf = e*sinh(Ff) - Ff;  % rad, "mean" anomaly (hyperbolic)
    % n = sqrt(mu/abs(a^3));              % rad/s, mean motion
    % tee = N_h0 / n;
    % F = funcF(dt, n, N_h0, e, xero);    % rad, *my* function used
end
% theta_t = 2*atan2( (sqrt(e+1)*tanh(F/2)), sqrt(e-1) ); % rad, new theta