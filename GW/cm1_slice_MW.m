%% Analysis of CM1 output: Obtain vertical slice (mountain waves)

% Output directory 
cm1_dir = './test_3D_hill';

% Variable name to take vertical slice
slice_var = {'uinterp', 'vinterp', 'winterp', 'thpert', 'prspert'};

%% Get variable names to slice vertical profiles

% Make directory to save data
slice_dir = fullfile(cm1_dir, 'vert_slice');
if exist(slice_dir, 'dir') ~= 7
    mkdir(slice_dir);
    disp(['Directory ', slice_dir, ' created successfully.']);
else
    disp(['Directory ', slice_dir, ' already exists.']);
end

% All NC output files
filelist = dir(fullfile(cm1_dir, 'cm1out_0*_i.nc'));
Nfile = length(filelist);

% Starting file number
i0 = sscanf(filelist(1).name, 'cm1out_%d_i.nc');

% Example file to obtain variable names & vertical grid
ref_file = fullfile(cm1_dir, ['cm1out_', sprintf('%06d', i0), '_i.nc']);
ref_info = ncinfo(ref_file);

% Coordinates
xh = ncread(ref_file, 'xh');  xf = ncread(ref_file, 'xf');
yh = ncread(ref_file, 'yh');  yf = ncread(ref_file, 'yf');
zh = ncread(ref_file, 'zh');  zf = ncread(ref_file, 'zf');

% Missing value
miss_value = ncreadatt(ref_file, "/", "missing_value");

% List of variable names (whose dimension is 4)
if any(strcmp(slice_var, 'all'))
    slice_var = {};
    for i = 1:length(ref_info.Variables)
        data_size = ref_info.Variables(i).Size;
    
        if length(data_size) == 4
            slice_var{end+1} = ref_info.Variables(i).Name;
        end
    end
end

% List of variable names
disp(['Number of variables: ', num2str(length(slice_var))]);
disp(strjoin(slice_var, ', '));

%% Vertical slice

% Create empty struct
slice_struct = struct;
time = zeros([1, Nfile]);

% Parallel setup
parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')));

% Vertical slices of variables
for i = 1:numel(slice_var) 
    varname = slice_var{i};
    nc_filepath = fullfile(cm1_dir, filelist(1).name);

    tmp = squeeze(ncread(nc_filepath, varname));
    sz = size(tmp);  iy_mid = ceil(sz(2)/4);
    var_mat = zeros([sz(1), sz(3), Nfile]);

    tmp = squeeze(tmp(:, iy_mid, :));
    tmp(tmp == miss_value) = NaN;
    var_mat(:, :, 1) = tmp;

    parfor j = 2:Nfile
        nc_filepath = fullfile(cm1_dir, filelist(j).name);
        tmp = squeeze(ncread(nc_filepath, varname));
        tmp = squeeze(tmp(:, iy_mid, :));
        tmp(tmp == miss_value) = NaN;
        var_mat(:, :, j) = tmp;
    end

    slice_struct.(varname) = var_mat;
end

% Time information
parfor j = 1:Nfile
    nc_filepath = fullfile(cm1_dir, filelist(j).name);
    time(j) = ncread(nc_filepath, 'time');
end
disp('Finish taking vertical slices.');

%% Save data

% Save data to matfile
y_slice = [yh(iy_mid), yf(iy_mid)];
save(fullfile(slice_dir, 'slice.mat'), "slice_struct", "time", ...
     "xh", "xf", "zh", "zf", "y_slice");
