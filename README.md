# GBM SPM12 L-TDI Pipeline

A MATLAB/SPM12 workflow for lesion-aware spatial normalization of glioblastoma MRI and label-wise calculation of the lesion–tract density index (L-TDI) using a normative streamline atlas.

## What this code does

1. Matches MRI images and multi-label lesion masks by filename.
2. Performs SPM12 Old Normalise with an inverse lesion weight image.
3. Resamples images and masks to a user-defined MNI grid.
4. Identifies normative atlas streamlines intersecting each requested lesion label.
5. Builds label-specific lesion–tract density maps (L-TDMs).
6. Computes L-TDI as the mean non-zero streamline density within each L-TDM.
7. Exports `ltdi_results.csv`, normalized images/masks, and optional L-TDM NIfTI files.

## Requirements

- MATLAB R2020b or later
- SPM12 with the Old Normalise toolbox available
- A 3D MNI template compatible with the source MRI modality
- A MAT-file containing `fibers_vox`, stored as a `1 × N` or `N × 1` cell array
- Paired 3D NIfTI images and masks

The code supports `.nii` and `.nii.gz` input files. Image and mask filenames must have identical stems:

```text
images/
├── case001.nii.gz
└── case002.nii.gz

masks/
├── case001.nii.gz
└── case002.nii.gz
```

## Installation

Clone the repository and open MATLAB in the repository directory:

```bash
git clone https://github.com/<YOUR-USERNAME>/gbm-spm12-ltdi.git
cd gbm-spm12-ltdi
```

Do not upload patient images, masks, local paths, or atlas files to a public repository.

## Configuration

Copy `config/gbm_spm12_ltdi_config_example.m` to a local file outside version control, then edit the paths and label definitions.

```matlab
addpath('src');
addpath('config');

cfg = gbm_spm12_ltdi_config_example();

cfg.spmDir = 'C:\path\to\spm12';
cfg.templatePath = 'C:\path\to\mni_template.nii';
cfg.imageDir = 'C:\path\to\images';
cfg.maskDir = 'C:\path\to\masks';
cfg.atlasMatPath = 'C:\path\to\dTOR_fibers_vox_1mm.mat';
cfg.outputDir = fullfile(pwd, 'outputs');

cfg.labels = [1 2 4];
cfg.labelNames = {'NecroticOrNonEnhancingCore', 'PeritumoralEdema', 'EnhancingTumor'};

resultsTable = gbm_spm12_ltdi_batch(cfg);
```

## Important configuration options

| Field | Meaning |
|---|---|
| `labels` | Mask values for which L-TDI will be calculated |
| `labelNames` | Human-readable names corresponding to `labels` |
| `runNormalization` | Set to `false` only when normalized files already exist in `outputs/normalized_images` and `outputs/normalized_masks` |
| `runLTDI` | Set to `false` to perform spatial normalization only |
| `saveLtdmNii` | Controls whether label-specific L-TDM NIfTI files are saved |
| `atlasReadMode` | `auto`, `load`, or `matfile`; use `matfile` for large atlases when memory is limited |
| `uniqueFiberVoxelsPerStreamline` | Removes repeated voxel visits within an individual streamline before density accumulation |

## Outputs

```text
outputs/
├── cfg_snapshot.mat
├── failed_cases.txt
├── ltdi_results.csv
├── normalized_images/
├── normalized_masks/
├── normalization_parameters/
└── ltdm_maps/
```

`ltdi_results.csv` includes case identifier, original and normalized file paths, label information, lesion voxel count, number of intersecting streamlines, non-zero L-TDM voxel count, and L-TDI.

## Data protection

This repository intentionally excludes raw imaging data, lesion masks, atlas files, intermediate files, and generated outputs. Before publication, confirm that your institutional data-sharing and software-licensing policies permit public release.

## Citation

If you use this code, please cite the accompanying article and the normative tract atlas used in your analysis. Add the final bibliographic reference and DOI to `CITATION.cff` before creating a release.

## License

This repository is distributed under the MIT License. Confirm that this license is compatible with your institutional and atlas-specific terms before public release.
