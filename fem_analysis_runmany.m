addpath("/d/gmi/1/simeon/simul_scalp/fieldtrip")
addpath("/d/gmi/1/simeon/simul_scalp/fieldtrip/external/spm12")
addpath("/d/gmi/1/simeon/hsc-eeg-tools/MATLAB")
ft_defaults

clear
load fieldtrip/template/headmodel/standard_seg.mat

cfg        = [];
cfg.shift  = 0.3;
cfg.method = 'hexahedral';
cfg.datatype = 'indexed';
hex_mesh = ft_prepare_mesh(cfg, mri);

% fill in missing labels
hex_mesh.tissuelabel = {'scalp', 'skull', 'brain'};

% load 10-20 electrodes
elec = ft_read_sens('standard_1020.elc');

% load sample pt localized electrodes
coords = readtable("/d/gmi/r1/CLAS/018/018.labels.csv");
coords = coords(strcmp(coords.Type, 'SEEG'), :);

% create original models
conductivities = struct('gray', 0.33, 'white', 0.14, 'csf', 1.79, 'brain', 1.5, 'skull', 0.01, 'scalp', 0.43, 'holes', 1.79);  % check brain

cfg        = [];
cfg.method = 'simbio';
cfg.conductivity = cellfun(@(x) conductivities.(x), hex_mesh.tissuelabel);

headmodel_orig  = ft_prepare_headmodel(cfg, hex_mesh);

cfg = [];
cfg.elec = elec;
cfg.resolution = 5;

cfg.headmodel = headmodel_orig;
sourcemodel_orig = ft_prepare_leadfield(cfg);
save('sourcemodel_orig.mat', 'sourcemodel_orig')


seeg_sens = [];
seeg_sens.label = coords.Label;
seeg_sens.chanpos = [coords.LocX, coords.LocY, coords.LocZ];
seeg_sens.elecpos = seeg_sens.chanpos;
seeg_sens.unit = 'mm';

% default parameters
% BURRHOLE_SIZE = 1.5;   % in mm
for BURRHOLE_SIZE = [0.5, 1.5, 3, 5, 10, 20]

    %% poke holes into meshes
    % get unique electrode names
    electrode_name = coords.Label(strcmp(coords.Type, 'SEEG'));
    for e = 1:length(electrode_name)
        tmp = regexp(electrode_name{e}, '^([A-Za-z]*)[0-9]*$', 'tokens');
        electrode_name{e} = tmp{1}{1};
    end
    unique_seeg = unique(electrode_name);
    
    hex_mesh_holey = hex_mesh;
    
    for uu = 1:length(unique_seeg)
        curr_chans = cellfun(@(x) startsWith(x, unique_seeg{uu}), coords.Label);
    
        % draw line through localized electrodes through to scalp
        line_interc = mean(seeg_sens.chanpos(curr_chans,:), 1);
        A = seeg_sens.chanpos(curr_chans,:) - line_interc;
        [U,S,~] = svd(A');
        line_slope = U(:,1);
        t = line_slope'*A';
        t1 = min(t);
        t2 = max(t);
    
        ln = line_interc' + [t1,t2] .* line_slope; % size 3x2
    
        % plot to check
    %     figure;
    %     hold on
    %     plot3(seeg_sens.chanpos(curr_chans,1), seeg_sens.chanpos(curr_chans,2), seeg_sens.chanpos(curr_chans,3), 'o')
    %     plot3(ln(1,:), ln(2,:), ln(3,:), 'r')
    % 
    %     title([unique_seeg{uu}, ' line fit'])
    
        % manually set where bone to CSF
        skull_idx = find(strcmp(hex_mesh.tissuelabel, 'skull'));
        if any(strcmp(hex_mesh.tissuelabel, 'csf'))
            csf_idx = find(strcmp(hex_mesh.tissuelabel, 'csf'));
        elseif any(strcmp(hex_mesh.tissuelabel, 'brain'))
            csf_idx = find(strcmp(hex_mesh.tissuelabel, 'brain'));
        else
            error('No brain-equivalent tissue type found');
        end
    
        % get middle of hex from position
        skh = find(hex_mesh.tissue == skull_idx);
    
        n_changed = 0;
        distances = zeros(size(hex_mesh.tissue));
    
        for kk = 1:length(skh)
            hex_pos = mean(hex_mesh.pos(hex_mesh.hex(skh(kk), :), :), 1);
    
            % how close to line
            b = line_interc - hex_pos;
            distance_to_line = norm(cross(line_slope, b)) / norm(line_slope);
            distances(skh(kk)) = distance_to_line;
    
            if distance_to_line < (BURRHOLE_SIZE / 2)
    %             hex_mesh_holey.tissue(skh(kk)) = csf_idx;
                hex_mesh_holey.tissue(skh(kk)) = 4;
                n_changed = n_changed + 1;
            end
        end
    
        fprintf('Changed %d hexes to csf.\n', n_changed)
        fprintf('Closest hex is %.3f mm away.\n', min(distances(distances > 0)))
    end
    
    hex_mesh_holey.tissuelabel{4} = 'holes';
    
    
    %% generate headmodels
    cfg        = [];
    cfg.method = 'simbio';
    cfg.conductivity = cellfun(@(x) conductivities.(x), hex_mesh_holey.tissuelabel);
    headmodel_hole  = ft_prepare_headmodel(cfg, hex_mesh_holey);


    %% generate leadfields
    cfg = [];
    cfg.elec = elec;
    cfg.resolution = 5;

    cfg.headmodel = headmodel_hole;
    sourcemodel_hole = ft_prepare_leadfield(cfg);
    save(sprintf('sourcemodel_hole_%.1f.mat', BURRHOLE_SIZE), 'sourcemodel_hole')

    dhole = cellfun(@(x) ~isempty(x), sourcemodel_hole.leadfield);
    dorig = cellfun(@(x) ~isempty(x), sourcemodel_orig.leadfield);
    
    dcommon = dhole & dorig;
    dcommon = find(dcommon);

    lf_orig = cat(3, sourcemodel_orig.leadfield{dcommon});
    lf_hole = cat(3, sourcemodel_hole.leadfield{dcommon});

    % compute leadfield differences
    rel_diff = zeros(length(dcommon), 1);
    for kk = 1:length(dcommon)
        lfo = sourcemodel_orig.leadfield{dcommon(kk)};
        lfh = sourcemodel_hole.leadfield{dcommon(kk)};
        rel_diff(kk) = mean(abs(lfo - lfh) ./ ((abs(lfo) + abs(lfh)) / 2), 'all');
    end

    cos_diff = zeros(length(dcommon), 1);
    for kk = 1:length(dcommon)
        lfo = sourcemodel_orig.leadfield{dcommon(kk)};
        lfh = sourcemodel_hole.leadfield{dcommon(kk)};
        cos_diff(kk) = abs(1 - (dot(lfo(:), lfh(:)) / (norm(lfo(:)) * norm(lfh(:)))));
    end

    norm_diff = zeros(length(dcommon), 1);
    for kk = 1:length(dcommon)
        lfo = sourcemodel_orig.leadfield{dcommon(kk)};
        lfh = sourcemodel_hole.leadfield{dcommon(kk)};
        norm_diff(kk) = rssq(lfo(:)) - rssq(lfh(:));
    end

    % compute distance to scalp
    dist_to_scalp = zeros(length(dcommon), 1);
    scalp_pts = hex_mesh.pos(hex_mesh.tissue == 1,:);
    for kk = 1:length(dcommon)
        dist_to_scalp(kk) = min(rssq(sourcemodel_orig.pos(dcommon(kk), :) - scalp_pts, 2));
    end

    save(sprintf('lfdata_%.1f.mat', BURRHOLE_SIZE), 'lf_orig', 'lf_hole', 'cos_diff', 'norm_diff', 'rel_diff', 'dist_to_scalp');

end

