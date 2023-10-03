addpath(fullfile(pwd, 'cbrewer2'));
addpath(fullfile(pwd, 'Robust_Statistical_Toolbox', 'graphic_functions'));
addpath(fullfile(pwd, 'RainCloudPlots', 'tutorial_matlab'));


load('sourcemodel_orig.mat', 'sourcemodel_orig')

%% main loop
BURRHOLE_SIZES = [0.5, 1.5, 3, 5, 10, 20];

grp_cos_diff = [];
grp_norm_diff = [];
for BURRHOLE_SIZE = BURRHOLE_SIZES
    load(sprintf('sourcemodel_hole_%.1f.mat', BURRHOLE_SIZE), 'sourcemodel_hole');

    dhole = cellfun(@(x) ~isempty(x), sourcemodel_hole.leadfield);
    dorig = cellfun(@(x) ~isempty(x), sourcemodel_orig.leadfield);

    dcommon = dhole & dorig;
    dcommon = find(dcommon);

    lf_orig = cat(3, sourcemodel_orig.leadfield{dcommon});
    lf_hole = cat(3, sourcemodel_hole.leadfield{dcommon});


    cos_diff = zeros(length(sourcemodel_orig.label), 1);
    for kk = 1:length(sourcemodel_orig.label)
        lfo = squeeze(rssq(lf_orig(kk,:,:),2));
        lfh = squeeze(rssq(lf_hole(kk,:,:),2));

        cos_diff(kk) = abs(1 - (dot(lfo(:), lfh(:)) / (norm(lfo(:)) * norm(lfh(:)))));
    end

    norm_diff = zeros(length(sourcemodel_orig.label), 1);
    for kk = 1:length(sourcemodel_orig.label)
        lfo = squeeze(rssq(lf_orig(kk,:,:),2));
        lfh = squeeze(rssq(lf_hole(kk,:,:),2));

        norm_diff(kk) = rssq(lfo(:)) - rssq(lfh(:));
    end



    % load cos diff into columns
    grp_cos_diff = [grp_cos_diff, cos_diff(:)];
    grp_norm_diff = [grp_norm_diff, norm_diff(:)];
end


%% compare different electrodes with intact skull
lf_orig = cat(3, sourcemodel_orig.leadfield{cellfun(@(x) ~isempty(x), sourcemodel_orig.leadfield)});

idx1 = find(strcmp(sourcemodel_orig.label, 'Fpz'));
idx2 = find(strcmp(sourcemodel_orig.label, 'Oz'));

lf1 = squeeze(rssq(lf_orig(idx1,:,:),2));
lf2 = squeeze(rssq(lf_orig(idx2,:,:),2));
cos_diff_sens = abs(1 - (dot(lf1(:), lf2(:)) / (norm(lf1(:)) * norm(lf2(:)))));
norm_diff_sens = abs(rssq(lf1(:)) - rssq(lf2(:)));

% preview
cos_diff_sens, norm_diff_sens



%% Make raincloud plots 
% box plot
hf = figure;
ax = subplot(1, 2, 1);
% boxplot(grp_cos_diff, Labels=arrayfun(@(x) sprintf('%.1f', x), BURRHOLE_SIZES, UniformOutput=false));
raincloud_plot(grp_cos_diff);
yline(cos_diff_sens)
ax.YLabel.String = 'Cosine difference';
ax.YLim = [0, max(max(grp_cos_diff(:)), cos_diff_sens) * 1.05];



ax = subplot(1, 2, 2);
% boxplot(grp_norm_diff, Labels=arrayfun(@(x) sprintf('%.1f', x), BURRHOLE_SIZES, UniformOutput=false));
yline(norm_diff_sens)
ax.YLabel.String = 'Norm difference';
ax.YLim = [0, max(max(grp_norm_diff(:)), norm_diff_sens) * 1.05];




