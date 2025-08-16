%% Analyze averaged vertical profile for several runs

% Get the screen size
screen = get(0, 'ScreenSize');
fig_pos = [0, 0, screen(3), screen(4)/2];

% Simulation to analyze 
sim_dir = './Isaac';
sim_list = {'Isaac_dx20m', 'Isaac_dx40m'};

% Legend for plots
sim_legend = {'dx = 20 m', 'dx = 40 m'};

% Time interval for recording diagnostic output
diagfrq = 60.0;

%% Universal variables

% Number of simulations
Ndir = length(sim_list);

% Create figure objects
fig1 = figure('Name', 'Wind', 'Position', fig_pos);  ax1 = cell(1,3);
for i = 1:length(ax1); ax1{i} = subplot(1,length(ax1),i); hold on; end
fig2 = figure('Name', 'Wind Variance', 'Position', fig_pos);  ax2 = cell(1,3);
for i = 1:length(ax2); ax2{i} = subplot(1,length(ax2),i); hold on; end
fig3 = figure('Name', 'theta & qv', 'Position', fig_pos);  ax3 = cell(1,2);
for i = 1:length(ax3); ax3{i} = subplot(1,length(ax3),i); hold on; end
fig4 = figure('Name', 'Vertical Flux', 'Position', fig_pos);  ax4 = cell(1,3);
for i = 1:length(ax4); ax4{i} = subplot(1,length(ax4),i); hold on; end
fig5 = figure('Name', 'Keff & Leff', 'Position', fig_pos);  ax5 = cell(1,2);
for i = 1:length(ax5); ax5{i} = subplot(1,length(ax5),i); hold on; end
fig6 = figure('Name', 'Velocity relations', 'Position', fig_pos);  ax6 = cell(1,2);
for i = 1:length(ax6); ax6{i} = subplot(1,length(ax6),i); hold on; end

%% Loop over simulations

for i = 1:length(sim_list)

    % Read profile
    load(fullfile(sim_dir, sim_list{i}, 'Data', 'diag_profile.mat'));
    load(fullfile(sim_dir, sim_list{i}, 'Data', 'evolution.mat'));

    % Time range for averaging
    istart = time_average(1) * 3600 / diagfrq + 1;
    iend = time_average(2) * 3600 / diagfrq + 1;

    % Output variables
    zh = scalar_struct.zh;  zf = wlev_struct.zf;
    u = scalar_struct.u;  v = scalar_struct.v;  wsp = scalar_struct.wsp;
    upup = wlev_struct.upup_w;  vpvp = wlev_struct.vpvp_w;  wpwp = wlev_struct.wpwp_w;
    ufr = wlev_struct.ufr + wlev_struct.ufd;  vfr = wlev_struct.vfr + wlev_struct.vfd;
    try
        ufs = wlev_struct.ufs + wlev_struct.ufw;  vfs = wlev_struct.vfs + wlev_struct.vfw;
    catch
        ufs = wlev_struct.ufs;  vfs = wlev_struct.vfs;
    end
    rtke = wlev_struct.rtke_w;  stke = wlev_struct.stke;
    rtau = sqrt(ufr.^2 + vfr.^2);  stau = sqrt(ufs.^2 + vfs.^2);
    ttau = sqrt((ufr + ufs).^2 + (vfr + vfs).^2);
    ttau_w = interp1(zf, ttau, zh, 'linear');

    % Effective eddy diffusivity & mixing length
    dudz = wlev_struct.dudz;  dvdz = wlev_struct.dvdz;
    rKeff = rtau ./ sqrt(dudz.^2 + dvdz.^2);
    try
        sKeff = wlev_struct.kmv + wlev_struct.kmw;
    catch
        sKeff = wlev_struct.kmv;
    end
    Leff = sqrt((rKeff+sKeff) ./ sqrt(dudz.^2 + dvdz.^2));

    % Obtain 10 m height variables
    dz = zh(2) - zh(1);  
    zh_mask = (zh <= 0.01 + dz*2/3) & (zh >= 0.01 - dz*2/3);
    zf_mask = (zf <= 0.01 + dz*2/3) & (zf >= 0.01 - dz*2/3);
    inflow_theta = -rad2deg(atan(mean(u(zh_mask) ./ v(zh_mask))));
    wsp10 = mean(wsp(zh_mask));

    % Normalization variables
    % ust = mean(evo_struct.ust(istart:iend));
    ust = sqrt(ttau(1));
    h_pbl = mean(evo_struct.hpbl(istart:iend));

    % Read input wind speed (hurr_vg)
    input_file = fullfile('./CM1_Output', sim_dir, 'Data', 'namelist.input');
    command = sprintf('grep "hurr_vg " %s | awk -F "=" ''{print $2}''', input_file);
    [~, result] = system(command);
    input_V = str2double(strtrim(result));
    if isnan(input_V)
        base_file = fullfile(sim_dir, sim_list{i}, 'Data', 'base_state.nc');
        base_v = squeeze(ncread(base_file, 'v'));
        input_V = base_v(2) - (base_v(2)-base_v(1))/dz * zh(2);
    end

    % Profile plots
    plot(ax1{1}, u, zh);  hold on;  plot(ax1{2}, v, zh);  hold on;  
    plot(ax1{3}, wsp, zh);  hold on;
    plot(ax2{1}, upup + 2/3.*stke, zf);  hold on;
    plot(ax2{2}, vpvp + 2/3.*stke, zf);  hold on;
    plot(ax2{3}, rtke + stke, zf);  hold on;
    plot(ax3{1}, scalar_struct.th, zh);  hold on;
    try
        obj = plot(ax3{2}, scalar_struct.qv*1e3, zh);  hold on;
        set(obj, 'DisplayName', sim_legend{i});
    catch
        warning('No moisture included in simulation: %s', sim_list{i});
    end
    plot(ax4{1}, ttau, zf);  hold on;  plot(ax4{2}, ttau./ust^2, zf);  hold on;
    plot(ax4{3}, ttau_w./wsp.^2, zh);  hold on;
    plot(ax5{1}, rKeff+sKeff, zf);  hold on;  plot(ax5{2}, Leff, zf);  hold on;
    
    % Scatter plots
    scatter(ax6{1}, input_V/wsp10, ust, 100, '+', 'LineWidth', 5);  hold on;
    scatter(ax6{2}, input_V, wsp10, 100, '+', 'LineWidth', 5);
end

%% Set figure properties

% Axes labels
xlabel(ax1{1}, '<u> (m/s)');                ylabel(ax1{1}, 'z (km)');
xlabel(ax1{2}, '<v> (m/s)');                ylabel(ax1{2}, 'z (km)');
xlabel(ax1{3}, '<U> (m/s)');                ylabel(ax1{3}, 'z (km)');
xlabel(ax2{1}, "<u'u'> (m^2/s^2)");         ylabel(ax2{1}, 'z (km)');
xlabel(ax2{2}, "<v'v'> (m^2/s^2)");         ylabel(ax2{2}, 'z (km)');
xlabel(ax2{3}, "TKE (m^2/s^2)");            ylabel(ax2{3}, 'z (km)');
xlabel(ax3{1}, "<\theta> (K)");             ylabel(ax3{1}, 'z (km)');
xlabel(ax3{2}, "<q_v> (g/kg)");             ylabel(ax3{2}, 'z (km)');
xlabel(ax4{1}, "|\tau|/\rho (m^2/s^2)");    ylabel(ax4{1}, 'z (km)');
xlabel(ax4{2}, "|\tau| / \rhou_*^2");       ylabel(ax4{2}, 'z (km)');
xlabel(ax4{3}, "|\tau| / \rho<U>^2");       ylabel(ax4{3}, 'z (km)');
xlabel(ax5{1}, "K_{eff} (m^2/s)");          ylabel(ax5{1}, 'z (km)');
xlabel(ax5{2}, "l_{eff} (m)");              ylabel(ax5{2}, 'z (km)');
xlabel(ax6{1}, "V / U10");                  ylabel(ax6{1}, 'u_* (m/s)');
xlabel(ax6{2}, "V (m/s)");                  ylabel(ax6{2}, 'U10 (m/s)');

% Axes limits
h_lim = [0, 3];

ylim(ax1{1}, h_lim);    ylim(ax1{2}, h_lim);    ylim(ax1{3}, h_lim);
ylim(ax2{1}, h_lim);    ylim(ax2{2}, h_lim);  
xlim(ax2{3}, [0, 40]);  ylim(ax2{3}, h_lim);
ylim(ax3{1}, [0, 2.5]); ylim(ax3{2}, [0, 2.5]);
ylim(ax4{1}, h_lim);
xlim(ax4{2}, [0, 1.1]); ylim(ax4{2}, h_lim);
ylim(ax4{3}, h_lim);     ylim(ax5{1}, h_lim);   ylim(ax5{2}, h_lim);

% Legends
legend(ax1{3}, sim_legend, 'Interpreter', 'none');
legend(ax2{3}, sim_legend, 'Interpreter', 'none');
legend(ax3{1}, sim_legend, 'Interpreter', 'none');
legend(ax3{2}, 'Interpreter', 'none');
legend(ax4{3}, sim_legend, 'Interpreter', 'none');
legend(ax5{2}, sim_legend, 'Interpreter', 'none');
legend(ax6{1}, sim_legend, 'Interpreter', 'none');
