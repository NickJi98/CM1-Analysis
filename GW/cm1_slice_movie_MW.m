%% Analysis of CM1 output: Movie of vertical slice (mountain waves)

addpath(genpath('/home/users/qingji/MATLAB'));

% Output directory 
cm1_dir = './test_3D_hill';

% Temporary figure directory
movie_dir = fullfile(cm1_dir, 'tmp_figs');

% Read vertical slice
slice_dir = fullfile(cm1_dir, 'vert_slice');
load(fullfile(slice_dir, 'slice.mat'));

% Variable names
slice_var = fieldnames(slice_struct);

% Number of frames
Nt = length(time);

% Colorbar: blue-white-red
mcolor = slanCM('bwr');

%% Movie of vertical slice

if exist(movie_dir, 'dir') ~= 7
    mkdir(movie_dir);
    disp(['Directory ', movie_dir, ' created successfully.']);
else
    disp(['Directory ', movie_dir, ' already exists.']);
end

Nxf = length(xf);  Nzf = length(zf);
for i = 1:numel(slice_var)
    
    % Vertical slice
    varname = slice_var{i};
    var_mat = slice_struct.(varname);

    % Coordinates
    sz = size(var_mat);
    xplot = xh;  if sz(1) == Nxf;  xplot = xf;  end
    zplot = zh;  if sz(2) == Nzf;  zplot = zf;  end

    % Color range
    cr = color_range(var_mat);

    % Plot vertical slice
    parfor j = 1:Nt
        figure('Visible', 'off');  colormap(mcolor);  hold on;
        set(gcf, 'Units', 'inches', 'Position', [1, 1, 7, 7]);
        pcolor(xplot, zplot, var_mat(:,:,j)');   shading interp;  clim(cr);
        xlabel('X (km)');   ylabel('Z (km)');
        xlim([xplot(1), xplot(end)]);  ylim([zplot(1), zplot(end)]);
        title(sprintf('%s, t = %d s', varname, time(j)), 'FontSize', 20);
        colorbar;

        ax = gca;
        exportgraphics(ax, fullfile(movie_dir, sprintf('%s_%04d.png', varname, j)), ...
            "Resolution", 200);
        close(gcf)
    end
end

%% Function: Determine color range

function cr = color_range(data)
    
    % Flattened array
    values = data(:);
    values = values(~isnan(values));
    
    % Basic statistics
    mean_val = mean(values);  rms_val  = std(values);
    
    % Check for symmetry about 0
    asymmetry = abs(mean_val) / rms_val;
    
    if asymmetry < 0.25
        % Symmetric about zero
        % bound = 0.9 * max(abs(values));
        bound = 2 * rms_val;
        cr = [-bound, bound];
    else
        % Centered around mean
        % bound = 0.9 * max(abs(values - mean_val));
        bound = 2 * rms_val;
        cr = [mean_val - bound, mean_val + bound];
    end
end
