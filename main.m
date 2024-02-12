% NovaSat Control calculations
mass = 20;%[kg]
dimensions = [20 20 30]/100; %[M]length width heights 
SAM = 500;                                  %[km] semi major axis
mu = 3.986e5;
orbit_altitude = 500                      ; %[km] 
Re = 6371                                 ; %[km] earths radius
R = Re+orbit_altitude;
latitude = 60*pi/180;                               %[rad]
Cd= 2.2 ;
cm = 0;% [M]center of mass
cp = 0.05;%center of pressure
cps = 0.05% solar center of pressure
Area = 0.2*0.3*2                          ; %[m^2]
Air_density = 10e-13;                       %[kg/m^3]
D = 0.5%[Am^2] residual dipole



lambda = asin( (SAM / R) * sin(latitude)) ; %magnetic latitude - ranges from 1 at equator to 2 at poles
B = 7.8e15*lambda/(R+orbit_altitude);       %[T]earths magnetic field
vel = sqrt(2*mu*(1/R - 1/(2*SAM)));%[km/s]
T = 2*pi*sqrt(SAM^3/mu);%[sec]

%% Inertia:
Izz = 1/12 * mass * (dimensions(2)^2+dimensions(1)^2);
Iyy = 1/12 * mass * (dimensions(1)^2+dimensions(3)^2);
Ixx = 1/12 * mass * (dimensions(2)^2+dimensions(3)^2);
%%Angular momentum and torque

%%atmospheric drag torque
Tad = 0.5*Air_density*Cd*Area*vel^2 * (cp-cm);
%%magnetic torque
Tm = D/B;
%%solar radiation toruqe
q = 0.6%solar coefficient
Tsr = 1367*Area*(1+q)*(cp-cm)
%%gravity gradient torque
GG_theta = 0%[rad]target of opprtunity angle
Tgg = 3*mu/(2*R^3) * abs(Izz-Iyy)*sin(2*GG_theta)

%%Torque from reaction wheels for distubance:
T_dominant = max([Tsr,Tad,Tm,Tgg])%[Nm]                           %enter the most dominant disturbance 
