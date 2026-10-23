# Software environment and validation scope

The source tree reflects the October 2026 GRIT workflows. Their complete production software versions
have not yet been recovered from GRIT. No March `renv.lock` or guessed version pin is presented as their
environment. Record `R.version`, package versions, Python packages, MATLAB version and cluster modules
from the actual run before claiming an environment reproduction.

## Static dependency inventory

R sources name `Matrix`, `R.matlab`, `arrow`, `data.table`, `dbscan`, `digest`, `doParallel`, `dplyr`,
`geosphere`, `glmnet`, `lubridate`, `matrixStats`, `mclust`, `pROC`, `qs`, `rlang`, `sensitivity`, `sf`,
`solitude`, `terra`, `testthat`, `tidyr`, `tidyselect`, `yaml` and `zoo`, in addition to base/recommended
packages. Some are optional or loaded indirectly. This is a source inventory, not a tested installation lock.

[requirements.txt](../requirements.txt) combines the Python imports and the post-processing dependency
list supplied with the scripts. It is intentionally unversioned until the production versions are established.
The OCIM calculation needs MATLAB and its external helper files; R parsing does not validate MATLAB execution.

Launchers also depend on Bash, Slurm commands, GNU command-line tools and site-specific module/library settings.
The checked-in partition/QOS, node exclusions, `R_LIBS_USER` and MATLAB executable defaults are GRIT examples
to review for another machine. `DFP_*` variables parameterise file locations, not the cluster scheduler.

## What was checked

Local source validation used R 4.5.3 and Python 3.12.14. The scientific payload contains 42 R, 16 Python and
32 shell/sbatch files, plus one MATLAB file. R was parsed without package loading; Python was compiled without
imports; shell files passed `bash -n`. The local translation identity check compared R/Python structures and
shell/MATLAB code after reversing explicitly declared transformations, and checked all 29 active MD5 guards.

The repository CI runs source-manifest checks, syntax checks and tests of the manifest checker. It does not
install scientific dependencies or run any AIS, spatial, sensitivity or OCIM calculation. Its selected R/Python
versions are syntax-validation tools, not assertions about the original GRIT environment.

An isolated execution trial on GRIT remains required: path existence and structural checks alone cannot prove
that a launcher never falls back to an older cluster copy or that every external input is available.
