%% Analysis of CM1 output: Make movies of vertical slice

% Output directory 
cm1_dir = './';

% Variable name to plot
movie_var = {'rhopert'};

% Color limit
crange = [-2e-5, 2e-5];

%% Get variable names to analyze vertical profiles

% Make directory to save data
movie_dir = fullfile(cm1_dir, 'tmp_figs');
if exist(movie_dir, 'dir') ~= 7
    mkdir(movie_dir);
    disp(['Directory ', movie_dir, ' created successfully.']);
else
    disp(['Directory ', movie_dir, ' already exists.']);
end

% All NC output files
filelist = dir(fullfile(cm1_dir, 'cm1out_0*.nc'));
Nfile = length(filelist);

% Starting file number
i0 = regexp(filelist(1).name, 'cm1out_(\d+).nc', 'tokens', 'once');
i0 = str2double(i0{1});

% Example file to obtain variable names & vertical grid
ref_file = fullfile(cm1_dir, ['cm1out_', sprintf('%06d', i0), '.nc']);
ref_info = ncinfo(ref_file);

% Coordinates
xh = ncread(ref_file, 'xh');  xf = ncread(ref_file, 'xf');
zh = ncread(ref_file, 'zh');  zf = ncread(ref_file, 'zf');
Nxh = length(xh);  Nxf = length(xf);
Nzh = length(zh);  Nzf = length(zf);

% List of variable names (whose dimension is 4)
if any(strcmp(movie_var, 'all'))
    movie_var = {};
    for i = 1:length(ref_info.Variables)
        data_size = ref_info.Variables(i).Size;
    
        if length(data_size) == 4
            movie_var{end+1} = ref_info.Variables(i).Name;
        end
    end
end

% List of variable names
disp(['Number of variables: ', num2str(length(movie_var))]);
disp(strjoin(movie_var, ', '));

%% Gravity wave dispersion

input_vars = {'var2', 'var9', 'var6', 'var19'};
vals = read_namelist(fullfile(cm1_dir, 'namelist.input'), input_vars);
xc = vals.var2 / 1e3;   zc = vals.var9 / 1e3;
omega = vals.var6;      N2 = vals.var19;
theta = asin(omega / sqrt(N2));

%% Plot vertical slice

% Parallel setup
parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')));

% Read variable
parfor j = 1:Nfile
    nc_filepath = fullfile(cm1_dir, filelist(j).name);
    iframe = regexp(filelist(j).name, 'cm1out_(\d+).nc', 'tokens', 'once');

    for i = 1:length(movie_var)
        varname = movie_var{i};
        var_mat = squeeze(ncread(nc_filepath, varname));

        sz = size(var_mat);
        xplot = xh;  if sz(1) == Nxf;  xplot = xf;  end
        zplot = zh;  if sz(3) == Nzf;  zplot = zf;  end
        var_mat = squeeze(var_mat(:, ceil(sz(2)/2), :));

        % Plot vertical slice
        figure('Visible', 'off');  colormap('gray');  hold on;
        set(gcf, 'Units', 'inches', 'Position', [1, 1, 7, 7]);
        pcolor(xplot, zplot, var_mat');   shading interp;  clim(crange);
        xlabel('X (km)');   ylabel('Y (km)');
        axis equal;  pbaspect([1,1,1]);
        title(sprintf('%s, frame %d', varname, str2double(iframe)), 'FontSize', 24);
        colorbar;

        % Group velocity direction
        plot(xh, zc + (xh-xc).*tan(theta), 'r-', 'LineWidth', 2);
        plot(xh, zc - (xh-xc).*tan(theta), 'r-', 'LineWidth', 2);
        xlim([xplot(1), xplot(end)]);  ylim([zplot(1), zplot(end)]);

        % Phase velocity direction
        alen = 0.5;
        quiver(1.5*xc, zc+0.5*xc*tan(theta), sin(theta)*alen, -cos(theta)*alen, 0, ...
            'b', 'LineWidth', 2, 'MaxHeadSize', 1.5);

        ax = gca;
        exportgraphics(ax, fullfile(movie_dir, sprintf('%s_%04d.png', varname, j)), ...
            "Resolution", 200);
    end
end

%% Function: Read input parameters

function extracted = read_namelist(filepath, vars_to_extract)
    extracted = struct();
    fid = fopen(filepath, 'r');
    if fid == -1
        error('Failed to open file: %s', filepath);
    end

    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if contains(line, '=')
            parts = strsplit(line, '=');
            var = strtrim(parts{1});
            val_str = strtrim(strtok(parts{2}, ','));

            if ismember(var, vars_to_extract)
                val_num = str2double(val_str);
                if isnan(val_num)
                    extracted.(var) = val_str;  % keep as string
                else
                    extracted.(var) = val_num;  % convert to number
                end
            end
        end
    end

    fclose(fid);
end
