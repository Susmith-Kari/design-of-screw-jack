clc; clear;

fprintf('**** Design of Screwjack ****\n\n');


W = input('Enter Load to be lifted (KN):');
Sst = input('Enter Screw material Tensile Strength (Mpa):');
Sss = input('Enter Screw material Shear Strength (Mpa):');
Pb = input('Enter Bearing pressure between nut and Screw (Mpa):');
Nf = input('Enter Factor of Safety:');
mu = input('Enter Coefficient of friction between screw and nut:');
mu1 = input('Enter Coefficient of collar friction:');
Snt = input('Enter Nut Tensile strength (Mpa):'); 
Snc = input('Enter Nut Compressive strength (Mpa):'); 
Sns = input('Enter Nut Shear strength (Mpa):');
max_lift = input('Enter Lift height (mm):');


fprintf('\n####----------Design of SCREW----------####\n');

% Safe allowable stress values
Sigma_S = Sst / Nf;  % Safe Tensile stress for screw
Tau_S = Sss / Nf;    % Safe Shear stress for screw

% Core Diameter Calculation
p = 10;  % Pitch of the screw
dc = sqrt((4 * W * 1000) / (pi * Sigma_S));  % Core diameter (dc)
do = round(1.6 * dc);   % Major diameter (do)
do1 = do + 2;
dc1 = do1 - p;
dm = (do1 + dc1) / 2;

% Displaying results
fprintf('Screw Major Diameter (do): %.2f mm\n', do1);
fprintf('Screw Minor Diameter (dc): %.2f mm\n', dc1);

% Helix and Friction Angle
q = p / (pi * dm);
alpha = atan(q);
phi = atan(mu);
alpha_deg = rad2deg(alpha);  % Helix angle in degrees
phi_deg = rad2deg(phi);      % Friction angle in degrees

fprintf('Screw Helix Angle: %.2f degrees\n', alpha_deg);
fprintf('Friction Angle: %.2f degrees\n', phi_deg);

% Torque required to rotate the screw
T1 = (W * 1000 * tan(alpha + phi) * (dm / 2));  % Torque in Nmm
fprintf('Torque required to rotate the screw T1: %.2f Nmm\n', T1);

% Torsional Shear stress due to torque T1
tau = (16 * T1) / (pi * (dc1^3));  
fprintf('Torsional Shear Stress: %.2f MPa\n', tau);

% Direct compressive stress due to axial load
sigmaC = (W * 1000) / (pi / 4 * dc1^2);
fprintf('Direct Compressive Stress: %.2f MPa\n', sigmaC);

% Safety check for shear and compressive stress
if tau < Tau_S && sigmaC < Sigma_S
    fprintf('Screw design is SAFE based on induced stresses.\n');
else
    fprintf('Screw design is NOT SAFE.\n');
end

% Maximum principal stresses
tau_max = sqrt(0.25 * (sigmaC^2 + 4 * tau^2));
sigma_max = 0.5 * sigmaC + tau_max;

fprintf('Max Shear Stress: %.2f MPa\n', tau_max);
fprintf('Max Normal Stress: %.2f MPa\n', sigma_max);

if tau_max < Tau_S && sigma_max < Sigma_S
    fprintf('Screw design is SAFE based on maximum principal stresses.\n');
else
    fprintf('Screw design is NOT SAFE.\n');
end

%----Design of Nut----%
fprintf('\n####----------Design of NUT----------####\n');

% Safe stresses for nut
Sigma_t_nut = Snt / Nf;  
Sigma_c_nut = Snc / Nf;  
Tau_nut_safe = Sns / Nf;

fprintf('Safe Tensile Stress for Nut: %.2f MPa\n', Sigma_t_nut);
fprintf('Safe Compressive Stress for Nut: %.2f MPa\n', Sigma_c_nut);
fprintf('Safe Shear Stress for Nut: %.2f MPa\n', Tau_nut_safe);

% Height of nut based on bearing pressure
bearing_area = (pi / 4) * (do1^2 - dc1^2);
n = ceil(W * 1000 / (Pb * bearing_area));  % Number of threads in contact
H = n * p;  % Height of the nut

fprintf('Height of Nut: %.2f mm\n', H);
fprintf('Number of threads in contact: %d\n', n);

% Transverse shear stress for screw and nut
t = p / 2;  % Thread thickness
Tau_screw = (W * 1000) / (pi * dc1 * n * t);
Tau_nut = (W * 1000) / (pi * do1 * n * t);

fprintf('Transverse Shear Stress for Screw: %.2f MPa\n', Tau_screw);
fprintf('Transverse Shear Stress for Nut: %.2f MPa\n', Tau_nut);

if Tau_screw < Tau_S && Tau_nut < Tau_nut_safe
    fprintf('Design of screw and nut is SAFE.\n');
else
    fprintf('Design of screw and nut is NOT SAFE.\n');
end

%----Nut Collar Design----%
fprintf('\n####----------Nut Collar Design----------####\n');

bb = (4 * W * 1000) / (pi * Sigma_t_nut);
D1 = ceil(sqrt(bb + do1^2));  % Collar inner diameter
fprintf('Collar Inner Diameter (D1): %.2f mm\n', D1);

vv = (4 * W * 1000) / (pi * Sigma_c_nut);
D2 = ceil(sqrt(vv + D1^2));  % Collar outer diameter
fprintf('Collar Outer Diameter (D2): %.2f mm\n', D2);

t1 = ceil(W * 1000 / (pi * D1 * Tau_nut_safe));  % Thickness of nut collar
fprintf('Nut Collar Thickness (t1): %.2f mm\n', t1);

%----Handle Design----%
fprintf('\n####----------Handle Dimensions----------####\n');
Rmean = (D3 / 2 + D4 / 2) / 2;
T2 = mu1 * W * 1000 * Rmean;  % Torque due to collar friction

T_handle = T1 + T2;  % Total torque
L_handle = T_handle / 300;  % Length of handle assuming 300 N force

fprintf('Length of Handle: %.2f m\n', L_handle / 1000);

Dh = ceil((32 * T_handle / (pi * 30))^(1/3));  % Diameter of handle
fprintf('Diameter of Handle (Dh): %.2f mm\n', Dh);

%----Body Design----%
fprintf('\n####----------Body Dimensions----------####\n');
D5 = 1.5 * D2;
t3 = 0.25 * do1;
D6 = ceil(2.25 * D2);
D7 = ceil(1.75 * D6);
Tb = 2 * t1;
Hb = ceil(max_lift + H + 1.75 * H);

fprintf('Body Diameter at Top (D5): %.2f mm\n', D5);
fprintf('Body Thickness (t3): %.2f mm\n', t3);
fprintf('Body Height (Hb): %.2f mm\n', Hb);

