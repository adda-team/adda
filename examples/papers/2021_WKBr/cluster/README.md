These scripts were used to produce the data for the paper using compute cluster and checkpoints to mitigate huge simulation times. It should be relatively easy to adapt them (in parts or as a whole) to other computational environments.

After setting the required parameters and paths to scripts (and uploading precomputed files with internal fields), execute
```
sh main_batch
```
on the cluster. This script will queue `batch` script several times to perform simulations with specific command line parameters.