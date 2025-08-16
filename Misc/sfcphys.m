%% Surface model parameters (Neutral Stability)

kappa = 0.4;    % Von Karman constant
z_ref = 10;     % Reference height (10 m)
S_ref = 14.88;   % Wind speed at reference height

%% Input: Roughness length

% Roughness length (m)
z0 = 0.29;

% Calculate other parameters
log_z = log(z_ref/z0 + 1);
Cd_ref = (kappa/log_z) ^ 2;
ust = S_ref * kappa / log_z;

fprintf('Input: z0 = %f m\n', z0);
fprintf('Cd = %f, ust = %f m/s\n\n', Cd_ref, ust);

%% Input: Friction velocity

% Friction velocity (m/s)
ust = 1.85;

% Calculate other parameters
log_z = S_ref * kappa / ust;
z0 = z_ref / (exp(log_z) - 1);
Cd_ref = (kappa/log_z) ^ 2;

fprintf('Input: ust = %f m/s\n', ust);
fprintf('Cd = %f, z0 = %f m\n\n', Cd_ref, z0);

%% Input: Drag coefficient

% Drag coefficient (at reference height)
Cd_ref = 0.013;

% Calculate other parameters
log_z = kappa / sqrt(Cd_ref);
z0 = z_ref / (exp(log_z) - 1);
ust = S_ref * kappa / log_z;

fprintf('Input: Cd = %f\n', Cd_ref);
fprintf('ust = %f m/s, z0 = %f m\n\n', ust, z0);
