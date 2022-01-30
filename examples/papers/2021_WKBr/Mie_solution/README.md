The script `exact-script.sh` calculates the exact internal electric field using either [`Bhfield`](https://scattport.org/index.php/light-scattering-software/mie-type-codes/list/521-bhfield) or [`Scattnlay`](https://github.com/ovidiopr/scattnlay) (to be obtained separately). The `toADDA.py` script changes Bhfield / Scattnlay grid so that it matches the ADDA one.

After setting the required parameters and paths to scripts, execute:
```
sh exact-script.sh
python toADDA.py
```
The scripts can also accept parameters through the command line - see the comments inside.