Exact_electric_field folder:
The script exact-script.sh calculates the exact internal electric field using either Bhfield or Scattnlay. You need to download the required project (Bhfield or Scattnlay) separately. In this case, the scripts find-lborder.py, find-rborder.py, find-radius.py are auxiliary for exact-script.sh. The toADDA.py script changes Bhfield / Scattnlay grid so that it matches the ADDA grid.

WKBr folder:
The scripts calculate WKBr internal electric field. For the calculation, you have to manually set the necessary parameters inside the main.py and run this script.

Cluster folder:
We used cluster to generate various data, including figures 15 and 16 from the paper. We used the scripts main_batch and batch. These scripts use checkpoints (read the ADDA manual). These scripts can be considered as an example for working with a cluster.

Otherwise, you can run Fig 15.sh and Fig 16.sh scripts that do not use checkpoints. It is assumed that computer / cluster resources for calculations are provided for an unlimited time.

All scripts are presented as examples and can help reproduce the results. In addition, in each case, you need to make sure that the scripts use the necessary parameters and paths to files. In other words, you need to manually configure the script before using it. However, despite this, we believe that scripts can be useful.
