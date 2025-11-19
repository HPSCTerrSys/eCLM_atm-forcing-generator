# Atmospheric forcing files based on other forcing

There exist a few global standard forcing data sets that can be downloaded together with their domain file from the ccsm data repository via this <a href="https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/atm/datm7/" target="_blank">link</a>. Based on this selection, it is easiest to start with these existing data files.

- <a href="https://escomp.github.io/ctsm-docs/versions/release-clm5.0/html/users_guide/setting-up-and-running-a-case/customizing-the-datm-namelist.html#clmgswp3v1-mode-and-it-s-datm-settings" target="_blank">GSWP3 NCEP forcing dataset</a> 
- <a href="https://escomp.github.io/ctsm-docs/versions/release-clm5.0/html/users_guide/setting-up-and-running-a-case/customizing-the-datm-namelist.html#clmcruncepv7-mode-and-it-s-datm-settings" target="_blank">CRUNCEP dataset</a> 
- <a href="https://escomp.github.io/ctsm-docs/versions/release-clm5.0/html/users_guide/setting-up-and-running-a-case/customizing-the-datm-namelist.html#clm-qian-mode-and-it-s-datm-settings" target="_blank">Qian dataset</a>

> [!CAUTION]
> The following part of the documentation is outdated!

For JSC users, an example python script to create forcings for a single-point case based on hourly observations can be found under `/p/scratch/cslts/shared_data/rlmod_eCLM`.

Simply copy the script to your directory and adapt it to your own data.

```sh
cp /p/scratch/cslts/shared_data/rlmod_eCLM/createnetCDF_forc_hourly_input.py $HOME
```
