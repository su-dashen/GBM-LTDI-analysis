# GitHub upload steps

1. Create a new empty repository on GitHub, for example `gbm-spm12-ltdi`.
2. Extract this package and open a terminal in the repository root.
3. Run the following commands:

```bash
git init
git add .
git commit -m "Initial public release"
git branch -M main
git remote add origin https://github.com/<YOUR-USERNAME>/gbm-spm12-ltdi.git
git push -u origin main
```

4. On GitHub, create a release tag matching the code used for the manuscript, for example `v1.0.0`.
5. Add the paper DOI and final article citation to `CITATION.cff` after publication.

Before pushing, run `git status` and confirm that no data, local configuration file, or output file is listed.
