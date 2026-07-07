function resultsTable = gbm_spm12_ltdi_batch(cfg)
if nargin < 1 || isempty(cfg)
    cfg = local_default_config();
else
    cfg = local_merge_struct(local_default_config(), cfg);
end

local_validate_config(cfg);
local_setup_spm(cfg);
local_ensure_dir(cfg.outputDir);
local_ensure_dir(cfg.normalizedImageDir);
local_ensure_dir(cfg.normalizedMaskDir);
local_ensure_dir(cfg.normalizationParamDir);
local_ensure_dir(cfg.tempDir);

if cfg.saveLtdmNii
    local_ensure_dir(cfg.ltdmDir);
end

pairs = local_collect_pairs(cfg);
if isempty(pairs)
    error('gbm_spm12_ltdi_batch:NoPairs', ...
        'No paired image and mask files were found.');
end

save(fullfile(cfg.outputDir, 'cfg_snapshot.mat'), 'cfg');

atlas = [];
if cfg.runLTDI
    atlas = local_prepare_atlas(cfg);
end

failedCasesFile = fullfile(cfg.outputDir, 'failed_cases.txt');
fidFailed = fopen(failedCasesFile, 'w');
if fidFailed == -1
    warning('gbm_spm12_ltdi_batch:FailedLog', ...
        'Unable to create the failed-case log: %s', failedCasesFile);
end
fileCloser = onCleanup(@() local_close_file(fidFailed));

caseRows = cell(height(pairs), 1);

for iCase = 1:height(pairs)
    caseId = pairs.case_id{iCase};
    fprintf('\n=== [%d/%d] %s ===\n', iCase, height(pairs), caseId);

    [normImageFile, normMaskFile, success] = local_prepare_case_mni(pairs(iCase, :), cfg);

    if ~success
        fprintf('Skipping case %s because normalisation failed.\n', caseId);
        if fidFailed ~= -1
            fprintf(fidFailed, '%s\n', caseId);
        end
        continue;
    end

    if ~cfg.runLTDI
        continue;
    end

    caseRows{iCase} = local_compute_case_ltdi( ...
        caseId, ...
        pairs.image_path{iCase}, ...
        pairs.mask_path{iCase}, ...
        normImageFile, ...
        normMaskFile, ...
        atlas, ...
        cfg);
end

validRows = caseRows(~cellfun(@isempty, caseRows));

if isempty(validRows)
    resultsTable = table();
else
    resultsTable = struct2table(vertcat(validRows{:}));
    writetable(resultsTable, fullfile(cfg.outputDir, 'ltdi_results.csv'));
end

clear fileCloser;

fprintf('\nFinished. Output directory: %s\n', cfg.outputDir);
fprintf('Failed cases: %s\n', failedCasesFile);
end


function cfg = local_default_config()
cfg = struct;

cfg.spmDir = '';
cfg.templatePath = '';
cfg.imageDir = '';
cfg.maskDir = '';
cfg.atlasMatPath = '';
cfg.outputDir = fullfile(pwd, 'outputs');

cfg.normalizedImageDir = fullfile(cfg.outputDir, 'normalized_images');
cfg.normalizedMaskDir = fullfile(cfg.outputDir, 'normalized_masks');
cfg.normalizationParamDir = fullfile(cfg.outputDir, 'normalization_parameters');
cfg.tempDir = fullfile(cfg.outputDir, 'temp');
cfg.ltdmDir = fullfile(cfg.outputDir, 'ltdm_maps');

cfg.labels = [1 2 4];
cfg.labelNames = {};
cfg.outputBoundingBox = [-90 -126 -72; 90 90 108];
cfg.outputVoxelSize = [1 1 1];
cfg.expectedMniDim = round( ...
    (cfg.outputBoundingBox(2, :) - cfg.outputBoundingBox(1, :)) ./ ...
    cfg.outputVoxelSize) + 1;

cfg.imageInterp = 4;
cfg.maskInterp = 0;
cfg.writePrefix = 'w';
cfg.overwrite = false;
cfg.cleanupTemp = false;
cfg.saveLtdmNii = true;

cfg.runNormalization = true;
cfg.runLTDI = true;

cfg.caseIds = {};
cfg.maxCases = inf;

cfg.atlasReadMode = 'auto';
cfg.atlasChunkSize = 20000;
cfg.maxFibers = inf;
cfg.uniqueFiberVoxelsPerStreamline = false;
cfg.progressEveryFibers = 100000;
end


function local_validate_config(cfg)
requiredFields = {'spmDir', 'imageDir', 'maskDir', 'outputDir'};

if cfg.runNormalization
    requiredFields{end + 1} = 'templatePath';
end

if cfg.runLTDI
    requiredFields{end + 1} = 'atlasMatPath';
end

for iField = 1:numel(requiredFields)
    fieldName = requiredFields{iField};
    if ~isfield(cfg, fieldName) || isempty(cfg.(fieldName))
        error('gbm_spm12_ltdi_batch:MissingConfig', ...
            'cfg.%s must be specified.', fieldName);
    end
end

if exist(cfg.spmDir, 'dir') ~= 7
    error('gbm_spm12_ltdi_batch:MissingSPM', ...
        'SPM directory not found: %s', cfg.spmDir);
end

if exist(cfg.imageDir, 'dir') ~= 7
    error('gbm_spm12_ltdi_batch:MissingImageDir', ...
        'Image directory not found: %s', cfg.imageDir);
end

if exist(cfg.maskDir, 'dir') ~= 7
    error('gbm_spm12_ltdi_batch:MissingMaskDir', ...
        'Mask directory not found: %s', cfg.maskDir);
end

if cfg.runNormalization && exist(cfg.templatePath, 'file') ~= 2
    error('gbm_spm12_ltdi_batch:MissingTemplate', ...
        'Template file not found: %s', cfg.templatePath);
end

if cfg.runLTDI && exist(cfg.atlasMatPath, 'file') ~= 2
    error('gbm_spm12_ltdi_batch:MissingAtlas', ...
        'Atlas file not found: %s', cfg.atlasMatPath);
end

if numel(cfg.labels) ~= numel(unique(cfg.labels))
    error('gbm_spm12_ltdi_batch:DuplicateLabels', ...
        'cfg.labels must contain unique label values.');
end

if ~isempty(cfg.labelNames) && numel(cfg.labelNames) ~= numel(cfg.labels)
    error('gbm_spm12_ltdi_batch:LabelNames', ...
        'cfg.labelNames must have the same length as cfg.labels.');
end
end


function local_setup_spm(cfg)
addpath(cfg.spmDir);
addpath(fullfile(cfg.spmDir, 'toolbox', 'OldNorm'));

if exist('spm', 'file') ~= 2 || exist('spm_normalise', 'file') ~= 2
    error('gbm_spm12_ltdi_batch:SPMSetup', ...
        'SPM or the OldNorm toolbox could not be located.');
end

spm('defaults', 'fmri');
spm_jobman('initcfg');
end


function pairs = local_collect_pairs(cfg)
imageFiles = local_list_nifti_files(cfg.imageDir);
maskFiles = local_list_nifti_files(cfg.maskDir);

imageMap = containers.Map('KeyType', 'char', 'ValueType', 'char');
maskMap = containers.Map('KeyType', 'char', 'ValueType', 'char');

for iFile = 1:numel(imageFiles)
    imageMap(local_strip_nifti_ext(imageFiles{iFile})) = imageFiles{iFile};
end

for iFile = 1:numel(maskFiles)
    maskMap(local_strip_nifti_ext(maskFiles{iFile})) = maskFiles{iFile};
end

caseIds = sort(intersect(keys(imageMap), keys(maskMap)));

if ~isempty(cfg.caseIds)
    requestedIds = cellstr(string(cfg.caseIds));
    caseIds = intersect(caseIds, requestedIds, 'stable');
end

if isfinite(cfg.maxCases)
    caseIds = caseIds(1:min(numel(caseIds), cfg.maxCases));
end

pairs = table( ...
    'Size', [numel(caseIds), 3], ...
    'VariableTypes', {'cell', 'cell', 'cell'}, ...
    'VariableNames', {'case_id', 'image_path', 'mask_path'});

for iCase = 1:numel(caseIds)
    caseId = caseIds{iCase};
    pairs.case_id{iCase} = caseId;
    pairs.image_path{iCase} = imageMap(caseId);
    pairs.mask_path{iCase} = maskMap(caseId);
end
end


function atlas = local_prepare_atlas(cfg)
atlas = struct;
atlas.readMode = lower(string(cfg.atlasReadMode));
atlas.chunkSize = cfg.atlasChunkSize;

switch atlas.readMode
    case "load"
        atlas.data = load(cfg.atlasMatPath, 'fibers_vox');
        [atlas.nFibers, atlas.orientation] = ...
            local_fiber_count_from_size(size(atlas.data.fibers_vox));
    case "matfile"
        atlas.data = matfile(cfg.atlasMatPath);
        [atlas.nFibers, atlas.orientation] = ...
            local_fiber_count_from_size(size(atlas.data, 'fibers_vox'));
    case "auto"
        try
            atlas.data = load(cfg.atlasMatPath, 'fibers_vox');
            atlas.readMode = "load";
            [atlas.nFibers, atlas.orientation] = ...
                local_fiber_count_from_size(size(atlas.data.fibers_vox));
            fprintf('Atlas read mode: load\n');
        catch loadError
            warning('gbm_spm12_ltdi_batch:AtlasFallback', ...
                ['Full atlas loading failed (%s). Falling back to matfile mode; ' ...
                'this may be slower.'], loadError.message);
            atlas.data = matfile(cfg.atlasMatPath);
            atlas.readMode = "matfile";
            [atlas.nFibers, atlas.orientation] = ...
                local_fiber_count_from_size(size(atlas.data, 'fibers_vox'));
            fprintf('Atlas read mode: matfile\n');
        end
    otherwise
        error('gbm_spm12_ltdi_batch:AtlasReadMode', ...
            'Unsupported cfg.atlasReadMode: %s', cfg.atlasReadMode);
end

atlas.nFibers = min(atlas.nFibers, cfg.maxFibers);
end


function [nFibers, orientation] = local_fiber_count_from_size(fiberSize)
if numel(fiberSize) < 2
    error('gbm_spm12_ltdi_batch:AtlasShape', ...
        'fibers_vox must be a vector-like cell array.');
end

if fiberSize(1) == 1
    nFibers = fiberSize(2);
    orientation = 'row';
elseif fiberSize(2) == 1
    nFibers = fiberSize(1);
    orientation = 'column';
else
    error('gbm_spm12_ltdi_batch:AtlasShape', ...
        'fibers_vox must be a 1-by-N or N-by-1 cell array.');
end
end


function [normImageFile, normMaskFile, success] = local_prepare_case_mni(caseRow, cfg)
normImageFile = '';
normMaskFile = '';
success = false;

caseId = caseRow.case_id{1};
normImageFile = fullfile(cfg.normalizedImageDir, sprintf('%s_mni.nii', caseId));
normMaskFile = fullfile(cfg.normalizedMaskDir, sprintf('%s_mni_labels.nii', caseId));
snMatFile = fullfile(cfg.normalizationParamDir, sprintf('%s_sn.mat', caseId));

if ~cfg.runNormalization
    if exist(normImageFile, 'file') ~= 2 || exist(normMaskFile, 'file') ~= 2
        error('gbm_spm12_ltdi_batch:MissingNormalizedFiles', ...
            ['cfg.runNormalization is false, but normalized image or mask files ' ...
            'are missing for case %s.'], caseId);
    end
    success = true;
    return;
end

if ~cfg.overwrite && ...
        exist(normImageFile, 'file') == 2 && ...
        exist(normMaskFile, 'file') == 2 && ...
        exist(snMatFile, 'file') == 2
    fprintf('Reusing normalized outputs for %s\n', caseId);
    success = true;
    return;
end

caseWorkDir = fullfile(cfg.tempDir, matlab.lang.makeValidName(caseId));
imageWorkDir = fullfile(caseWorkDir, 'image');
maskWorkDir = fullfile(caseWorkDir, 'mask');

local_ensure_dir(caseWorkDir);
local_ensure_dir(imageWorkDir);
local_ensure_dir(maskWorkDir);

sourceImageNii = local_materialize_nifti( ...
    caseRow.image_path{1}, imageWorkDir, 'source_image', cfg.overwrite);
sourceMaskNii = local_materialize_nifti( ...
    caseRow.mask_path{1}, maskWorkDir, 'source_mask', cfg.overwrite);
sourceWeightNii = fullfile(caseWorkDir, 'source_weight.nii');

try
    local_validate_pair_geometry(sourceImageNii, sourceMaskNii);
    local_write_inverse_weight_mask(sourceMaskNii, sourceWeightNii, cfg.overwrite);

    estimateFlags = spm_get_defaults('old.normalise.estimate');
    estimateFlags.graphics = 0;

    templateVol = spm_vol(cfg.templatePath);
    sourceVol = spm_vol(sourceImageNii);

    spm_normalise( ...
        templateVol, ...
        sourceVol, ...
        snMatFile, ...
        '', ...
        sourceWeightNii, ...
        estimateFlags);

    spm_write_sn( ...
        spm_vol(sourceImageNii), ...
        snMatFile, ...
        local_get_write_flags(cfg, cfg.imageInterp));

    spm_write_sn( ...
        spm_vol(sourceMaskNii), ...
        snMatFile, ...
        local_get_write_flags(cfg, cfg.maskInterp));

    tmpNormImage = fullfile( ...
        imageWorkDir, sprintf('%s%s.nii', cfg.writePrefix, 'source_image'));
    tmpNormMask = fullfile( ...
        maskWorkDir, sprintf('%s%s.nii', cfg.writePrefix, 'source_mask'));

    if exist(tmpNormImage, 'file') ~= 2 || exist(tmpNormMask, 'file') ~= 2
        error('gbm_spm12_ltdi_batch:MissingSPMOutput', ...
            'SPM did not generate both normalized output files for case %s.', caseId);
    end

    local_safe_copy(tmpNormImage, normImageFile, cfg.overwrite);
    local_safe_copy(tmpNormMask, normMaskFile, cfg.overwrite);
    local_validate_normalized_mask(normMaskFile, cfg);

    success = true;
catch caseError
    fprintf('Normalization failed for %s: %s\n', caseId, caseError.message);
end

if cfg.cleanupTemp && exist(sourceWeightNii, 'file') == 2
    delete(sourceWeightNii);
end
end


function rows = local_compute_case_ltdi( ...
    caseId, imageFile, maskFile, normImageFile, normMaskFile, atlas, cfg)

Vmask = spm_vol(normMaskFile);
maskData = round(spm_read_vols(Vmask));
maskData(~isfinite(maskData)) = 0;

if ~isequal(Vmask.dim(1:3), cfg.expectedMniDim)
    error('gbm_spm12_ltdi_batch:MaskGrid', ...
        ['Normalized mask dimensions do not match the expected MNI grid ' ...
        'for case %s.'], caseId);
end

nLabels = numel(cfg.labels);
densityMaps = cell(1, nLabels);
lesionVoxelCounts = zeros(1, nLabels);
intersectingFibers = zeros(1, nLabels);

for iLabel = 1:nLabels
    densityMaps{iLabel} = zeros( ...
        Vmask.dim(1), Vmask.dim(2), Vmask.dim(3), 'uint32');
    lesionVoxelCounts(iLabel) = nnz(maskData == cfg.labels(iLabel));
end

labelUnion = maskData > 0;

fprintf('Computing L-TDI for %s across %d fibers...\n', caseId, atlas.nFibers);

for firstFiber = 1:atlas.chunkSize:atlas.nFibers
    lastFiber = min(firstFiber + atlas.chunkSize - 1, atlas.nFibers);
    fiberChunk = local_get_atlas_chunk(atlas, firstFiber, lastFiber);

    for iFiber = 1:numel(fiberChunk)
        coords = fiberChunk{iFiber};

        if isempty(coords)
            continue;
        end

        coords = local_filter_valid_voxels( ...
            round(coords), Vmask.dim(1:3));

        if isempty(coords)
            continue;
        end

        lin = sub2ind( ...
            Vmask.dim(1:3), ...
            coords(:, 1), ...
            coords(:, 2), ...
            coords(:, 3));

        if cfg.uniqueFiberVoxelsPerStreamline
            lin = unique(lin, 'stable');
        end

        if ~any(labelUnion(lin))
            continue;
        end

        hitLabels = unique(maskData(lin));
        hitLabels(hitLabels == 0) = [];
        [isRequested, requestedIndices] = ismember(hitLabels, cfg.labels);
        requestedIndices = requestedIndices(isRequested);

        for iHit = 1:numel(requestedIndices)
            labelIndex = requestedIndices(iHit);
            densityMaps{labelIndex}(lin) = densityMaps{labelIndex}(lin) + 1;
            intersectingFibers(labelIndex) = ...
                intersectingFibers(labelIndex) + 1;
        end
    end

    if lastFiber == atlas.nFibers || ...
            mod(lastFiber, cfg.progressEveryFibers) == 0
        fprintf('[%s] fibers processed: %d / %d\n', ...
            caseId, lastFiber, atlas.nFibers);
    end
end

rows = repmat(struct( ...
    'case_id', '', ...
    'image_file', '', ...
    'mask_file', '', ...
    'normalized_image_file', '', ...
    'normalized_mask_file', '', ...
    'label', 0, ...
    'label_name', '', ...
    'lesion_voxels', 0, ...
    'intersecting_fibers', 0, ...
    'ltdm_nonzero_voxels', 0, ...
    'ltdi', NaN), nLabels, 1);

for iLabel = 1:nLabels
    densityMap = densityMaps{iLabel};
    nonzeroMask = densityMap > 0;
    nonzeroVoxelCount = nnz(nonzeroMask);

    if nonzeroVoxelCount == 0
        ltdi = NaN;
    else
        ltdi = mean(double(densityMap(nonzeroMask)));
    end

    if cfg.saveLtdmNii
        ltdmFile = fullfile( ...
            cfg.ltdmDir, ...
            sprintf('%s_label%d_ltdm.nii', caseId, cfg.labels(iLabel)));
        local_write_volume_like(Vmask, ltdmFile, single(densityMap), 'float32');
    end

    rows(iLabel).case_id = caseId;
    rows(iLabel).image_file = imageFile;
    rows(iLabel).mask_file = maskFile;
    rows(iLabel).normalized_image_file = normImageFile;
    rows(iLabel).normalized_mask_file = normMaskFile;
    rows(iLabel).label = cfg.labels(iLabel);
    rows(iLabel).label_name = local_get_label_name(cfg, iLabel);
    rows(iLabel).lesion_voxels = lesionVoxelCounts(iLabel);
    rows(iLabel).intersecting_fibers = intersectingFibers(iLabel);
    rows(iLabel).ltdm_nonzero_voxels = nonzeroVoxelCount;
    rows(iLabel).ltdi = ltdi;
end
end


function labelName = local_get_label_name(cfg, labelIndex)
if isempty(cfg.labelNames)
    labelName = sprintf('label%d', cfg.labels(labelIndex));
else
    labelName = char(string(cfg.labelNames{labelIndex}));
end
end


function writeFlags = local_get_write_flags(cfg, interpolation)
writeFlags = spm_get_defaults('old.normalise.write');
writeFlags.preserve = 0;
writeFlags.bb = cfg.outputBoundingBox;
writeFlags.vox = cfg.outputVoxelSize;
writeFlags.interp = interpolation;
writeFlags.wrap = [0 0 0];
writeFlags.prefix = cfg.writePrefix;
end


function local_validate_pair_geometry(imageNii, maskNii)
Vimage = spm_vol(imageNii);
Vmask = spm_vol(maskNii);

if numel(Vimage) ~= 1 || numel(Vmask) ~= 1
    error('gbm_spm12_ltdi_batch:Dimensionality', ...
        'Only 3D NIfTI files are supported.');
end

if ~isequal(Vimage.dim(1:3), Vmask.dim(1:3))
    error('gbm_spm12_ltdi_batch:ImageMaskDimensions', ...
        'Image and mask dimensions do not match.');
end

if max(abs(Vimage.mat(:) - Vmask.mat(:))) > 1e-4
    error('gbm_spm12_ltdi_batch:ImageMaskAffine', ...
        'Image and mask affine matrices do not match.');
end
end


function local_validate_normalized_mask(normMaskFile, cfg)
Vmask = spm_vol(normMaskFile);

if ~isequal(Vmask.dim(1:3), cfg.expectedMniDim)
    error('gbm_spm12_ltdi_batch:NormalizedMaskGrid', ...
        'Normalized mask does not match the expected MNI grid: %s', normMaskFile);
end
end


function local_write_inverse_weight_mask(maskNii, weightNii, overwrite)
if ~overwrite && exist(weightNii, 'file') == 2
    return;
end

Vmask = spm_vol(maskNii);
maskData = spm_read_vols(Vmask);
maskData(~isfinite(maskData)) = 0;
weightData = single(maskData == 0);

local_write_volume_like(Vmask, weightNii, weightData, 'float32');
end


function local_write_volume_like(Vref, outFile, data, dtypeName)
Vout = Vref;
Vout.fname = outFile;
Vout.dt = [spm_type(dtypeName) spm_platform('bigend')];
Vout.pinfo = [1; 0; 0];

spm_write_vol(Vout, data);
end


function outFile = local_materialize_nifti(inFile, outDir, outStem, overwrite)
[~, ~, ext] = fileparts(inFile);

if strcmpi(ext, '.gz')
    outFile = fullfile(outDir, sprintf('%s.nii', outStem));
else
    outFile = fullfile(outDir, sprintf('%s%s', outStem, ext));
end

if ~overwrite && exist(outFile, 'file') == 2
    return;
end

if strcmpi(ext, '.gz')
    unpackedFiles = gunzip(inFile, outDir);
    movefile(unpackedFiles{1}, outFile, 'f');
else
    copyfile(inFile, outFile, 'f');
end
end


function coords = local_filter_valid_voxels(coords, volumeSize)
isValid = ...
    coords(:, 1) >= 1 & coords(:, 1) <= volumeSize(1) & ...
    coords(:, 2) >= 1 & coords(:, 2) <= volumeSize(2) & ...
    coords(:, 3) >= 1 & coords(:, 3) <= volumeSize(3);

coords = coords(isValid, :);
end


function fiberChunk = local_get_atlas_chunk(atlas, firstFiber, lastFiber)
switch atlas.readMode
    case "load"
        if strcmp(atlas.orientation, 'row')
            fiberChunk = atlas.data.fibers_vox(1, firstFiber:lastFiber);
        else
            fiberChunk = atlas.data.fibers_vox(firstFiber:lastFiber, 1);
        end
    case "matfile"
        if strcmp(atlas.orientation, 'row')
            fiberChunk = atlas.data.fibers_vox(1, firstFiber:lastFiber);
        else
            fiberChunk = atlas.data.fibers_vox(firstFiber:lastFiber, 1);
        end
    otherwise
        error('gbm_spm12_ltdi_batch:AtlasReadMode', ...
            'Unsupported atlas read mode: %s', atlas.readMode);
end
end


function files = local_list_nifti_files(rootDir)
niiFiles = dir(fullfile(rootDir, '*.nii'));
niiGzFiles = dir(fullfile(rootDir, '*.nii.gz'));
fileInfo = [niiFiles; niiGzFiles];

files = cell(numel(fileInfo), 1);
for iFile = 1:numel(fileInfo)
    files{iFile} = fullfile(fileInfo(iFile).folder, fileInfo(iFile).name);
end
end


function stem = local_strip_nifti_ext(filePath)
[~, name, ext] = fileparts(filePath);

if strcmpi(ext, '.gz')
    [~, stem, ~] = fileparts(name);
else
    stem = name;
end
end


function local_safe_copy(sourceFile, destinationFile, overwrite)
if overwrite && exist(destinationFile, 'file') == 2
    delete(destinationFile);
end

copyfile(sourceFile, destinationFile, 'f');
end


function local_ensure_dir(dirPath)
if exist(dirPath, 'dir') ~= 7
    mkdir(dirPath);
end
end


function local_close_file(fileId)
if fileId ~= -1
    fclose(fileId);
end
end


function out = local_merge_struct(base, override)
out = base;

if isempty(override)
    return;
end

overrideFields = fieldnames(override);
for iField = 1:numel(overrideFields)
    fieldName = overrideFields{iField};
    out.(fieldName) = override.(fieldName);
end

if isfield(override, 'outputDir')
    if ~isfield(override, 'normalizedImageDir')
        out.normalizedImageDir = fullfile(out.outputDir, 'normalized_images');
    end

    if ~isfield(override, 'normalizedMaskDir')
        out.normalizedMaskDir = fullfile(out.outputDir, 'normalized_masks');
    end

    if ~isfield(override, 'normalizationParamDir')
        out.normalizationParamDir = ...
            fullfile(out.outputDir, 'normalization_parameters');
    end

    if ~isfield(override, 'tempDir')
        out.tempDir = fullfile(out.outputDir, 'temp');
    end

    if ~isfield(override, 'ltdmDir')
        out.ltdmDir = fullfile(out.outputDir, 'ltdm_maps');
    end
end

if isfield(override, 'outputBoundingBox') || ...
        isfield(override, 'outputVoxelSize')
    out.expectedMniDim = round( ...
        (out.outputBoundingBox(2, :) - out.outputBoundingBox(1, :)) ./ ...
        out.outputVoxelSize) + 1;
end
end
