projectRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(projectRoot, 'src'));
addpath(fullfile(projectRoot, 'config'));

cfg = gbm_spm12_ltdi_config_example();

cfg.spmDir = 'C:\path\to\spm12';
cfg.templatePath = 'C:\path\to\mni_template.nii';
cfg.imageDir = 'C:\path\to\images';
cfg.maskDir = 'C:\path\to\masks';
cfg.atlasMatPath = 'C:\path\to\dTOR_fibers_vox_1mm.mat';
cfg.outputDir = fullfile(projectRoot, 'outputs');

resultsTable = gbm_spm12_ltdi_batch(cfg);
disp(resultsTable);
