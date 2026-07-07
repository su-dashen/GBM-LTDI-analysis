function cfg = gbm_spm12_ltdi_config_example()
cfg = struct;

cfg.spmDir = 'C:\path\to\spm12';
cfg.templatePath = 'C:\path\to\mni_template.nii';
cfg.imageDir = 'C:\path\to\images';
cfg.maskDir = 'C:\path\to\masks';
cfg.atlasMatPath = 'C:\path\to\dTOR_fibers_vox_1mm.mat';
cfg.outputDir = fullfile(pwd, 'outputs');

cfg.labels = [1 2 4];
cfg.labelNames = {'Label1', 'Label2', 'Label4'};

cfg.outputBoundingBox = [-90 -126 -72; 90 90 108];
cfg.outputVoxelSize = [1 1 1];

cfg.imageInterp = 4;
cfg.maskInterp = 0;
cfg.writePrefix = 'w';

cfg.runNormalization = true;
cfg.runLTDI = true;
cfg.saveLtdmNii = true;
cfg.overwrite = false;
cfg.cleanupTemp = false;

cfg.atlasReadMode = 'auto';
cfg.atlasChunkSize = 20000;
cfg.maxFibers = inf;
cfg.uniqueFiberVoxelsPerStreamline = false;
cfg.progressEveryFibers = 100000;

cfg.caseIds = {};
cfg.maxCases = inf;
end
