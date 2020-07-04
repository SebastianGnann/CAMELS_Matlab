# CAMELS_Matlab
Loads CAMELS, CAMELS-GB, and CAMELS-BR data from original files and saves data as Matlab struc files.
CAMELS, CAMELS-CL, CAMELS-GB (TO DO: update to new format), and CAMELS-BR data are saved locally.

To limit file size, only one P and PET product is loaded:
- for CAMELS US, PET is calculated using a calibrated Priestley-Taylor coefficient which is now set to the default value of 1.26
- for CAMELS CL, MSWEP P and HARGREAVES PET are used
- for CAMELS GB, PET without interception correction is used
- for CAMELS BR, CHIRPS P, GLEAM PET, and CPC T are used

## Papers describing datasets:

Addor, N., Newman, A.J., Mizukami, N. and Clark, M.P., 2017. The CAMELS data set: catchment attributes and meteorology for large-sample studies. Hydrology and Earth System Sciences (HESS), 21(10), pp.5293-5313.

Alvarez-Garreton, C., Mendoza, P.A., Boisier, J.P., Addor, N., Galleguillos, M., Zambrano-Bigiarini, M., Lara, A., Cortes, G., Garreaud, R., McPhee, J. and Ayala, A., 2018. The CAMELS-CL dataset: catchment attributes and meteorology for large sample studies-Chile dataset. Hydrology and Earth System Sciences, 22(11), pp.5817-5846.

Coxon, G., Addor, N., Bloomfield, J.P., Freer, J., Fry, M. and Hannaford, J., 2020. CAMELS-GB: Hydrometeorological time series and landscape attributes for 671 catchments in Great Britain. Earth Syst. Sci. Data Discuss, 2020, pp.1-34.

Chagas, V.B., Chaffe, P.L., Addor, N., Fan, F.M., Fleischmann, A.S., Paiva, R.C. and Siqueira, V.A., 2020. CAMELS-BR: Hydrometeorological time series and landscape attributes for 897 catchments in Brazil. Earth System Science Data Discussions, pp.1-41.


## Links to datasets:

https://ral.ucar.edu/solutions/products/camels

https://doi.pangaea.de/10.1594/PANGAEA.894885?format=html#download

https://catalogue.ceh.ac.uk/documents/8344e4f3-d2ea-44f5-8afa-86d2987543a9

https://zenodo.org/record/3709338#.Xv4ciShKhPY


