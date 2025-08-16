%% Space-time correlation function

% Get the screen size
screen = get(0, 'ScreenSize');

% Colorbar: red-white-blue
mcolor = slanCM('bwr');  mcolor = mcolor(end-16:-1:16, :);

%% Load data

%%% Data file %%%
cm1_dir = './Isaac/R105km_dx20m_h2';
cm1_file = fullfile(cm1_dir, 'cm1out.nc');
fprintf('Read in data: %s\n\n', cm1_file);

%%% Read pressure source %%%
src = read_cm1_output(cm1_file);

% Mesh size (km, s)
dx = src.xh(2) - src.xh(1);  dy = src.yh(2) - src.yh(1);
dt = src.time(2) - src.time(1);

% Data format: Nx * Ny * Nt = 256 * 256 * 3600
[Nx, Ny, Nt] = size(src.pp);

%% Rotation of axis

% Rotation angle (anticlockwise in degree)
theta = 116;
fname = sprintf('corr_xt_%ddeg.mat', theta);

% Initialize array
pp_rot = zeros(Nx,Ny,Nt);

% Periodic extension of 2D field: Parameters
sx = Nx/4;  sy = Ny/4;
x_ex = (1:Nx+sx*2).*dx;     y_ex = (1:Ny+sy*2).*dy;
x_ex = x_ex - mean(x_ex);   y_ex = y_ex - mean(y_ex);

for it = 1:Nt

    % Periodic extension of 2D field
    pp = src.pp(:,:,it);
    pp_ex = [pp(end-sx+1:end, end-sy+1:end), pp(end-sx+1:end, :), pp(end-sx+1:end, 1:sy);
             pp(:, end-sy+1:end), pp, pp(:, 1:sy);
             pp(1:sx, end-sy+1:end), pp(1:sx, :), pp(1:sx, 1:sy)];

    % Interpolant
    [X, Y] = ndgrid(x_ex, y_ex);
    F_interp = griddedInterpolant(X, Y, pp_ex, 'linear', 'none');
    clear X Y pp_ex;

    % Rotated axis
    [X, Y] = ndgrid(src.xh-mean(src.xh), src.yh-mean(src.yh));
    X_rot = X .* cosd(theta) - Y .* sind(theta);
    Y_rot = X .* sind(theta) + Y .* cosd(theta);
    pp_rot(:,:,it) = F_interp(X_rot, Y_rot);
    clear X Y X_rot Y_rot;

end

%% Space-time correlations

% Max offset index in x and t
off_x_max = 32;  off_t_max = 32;

% Max domain for calculating correlation
calc_x_max = 128;  calc_t_max = 128;

% Lag space and time axes (km, s)
x_lag = (0:1:off_x_max).*dx;  t_lag = (0:1:off_t_max).*dt;

% Correlation matrix
Nx_ = off_x_max+1;  Ny_ = Nx_;  Nt_ = off_t_max+1;
corr_xt = zeros(Nx_, Nt_);
corr_yt = zeros(Ny_, Nt_);

% Calculate correlation
for it = 1:Nt_
    disp(it);
    for ix = 1:Nx_
        corr_xt(ix,it) = mean(pp_rot(1:calc_x_max, :, 1:calc_t_max) .* ...
            pp_rot(ix:ix+calc_x_max-1, :, it:it+calc_t_max-1), 'all');
    end
end

% Calculate correlation
for it = 1:Nt_
    for iy = 1:Ny_
        corr_yt(iy,it) = mean(pp_rot(:, 1:calc_x_max, 1:calc_t_max) .* ...
            pp_rot(:, iy:iy+calc_x_max-1, it:it+calc_t_max-1), 'all');
    end
end

% Save to mat file
save(fname, 'corr_xt', 'corr_yt', 'x_lag', 't_lag');

%% Convection velocity (x-direction)

load(fname);

% Ridge of maximum
corr_xt = corr_xt ./ max(corr_xt, [], 'all');
[corr_max, it_max] = max(corr_xt, [], 2);  t_max = t_lag(it_max);

% Linear fit
mask = (corr_max >= exp(-1));
P = polyfit(x_lag(mask), t_max(mask), 1);
Uc = 1/P(1);
fprintf('Convective velocity (%d deg): %g m/s\n', theta, Uc*1e3);

% Plot space-time correlation
figure('Name', 'Correlation');
colormap('jet');
contourf(x_lag, t_lag, corr_xt');  shading interp;  hold on;
plot(x_lag, P(1).*x_lag+P(2), 'w-', 'LineWidth', 3);
scatter(x_lag, t_max, 20, 'k', 'filled');
xlabel('$$r$$ (km)', 'interpreter', 'latex');   xlim([0, 0.5]);
ylabel('$$\tau$$ (s)', 'interpreter', 'latex'); ylim([0, 30]);
pbaspect([1.5,1,1]);  colorbar;

%% Convection velocity (y-direction)

load(fname);

% Ridge of maximum
corr_xt = corr_xt ./ max(corr_xt, [], 'all');
[~, it_max] = max(corr_yt, [], 2);  t_max = t_lag(it_max);

% Linear fit
mask = (x_lag < 0.51);
P = polyfit(x_lag(mask), t_max(mask), 1);
Uc = 1/P(1);
fprintf('Convective velocity (%d deg): %g m/s\n', theta-90, Uc*1e3);

% Plot space-time correlation
figure('Name', 'Correlation');
colormap('jet');
contourf(x_lag, t_lag, corr_yt');  shading interp; hold on;
plot(x_lag, P(1).*x_lag+P(2), 'w-', 'LineWidth', 3);
scatter(x_lag, t_max, 20, 'k', 'filled');
xlabel('$$r$$ (km)', 'interpreter', 'latex');   xlim([0, 0.5]);
ylabel('$$\tau$$ (s)', 'interpreter', 'latex'); ylim([0, 30]);
pbaspect([1.5,1,1]);  colorbar;
