# CAMELS_Matlab
Loads CAMELS, CAMELS-CL, CAMELS-GB, CAMELS-BR, and CAMELS-AUS data and saves the data as Matlab struct files.
The data and the resulting struct files are saved locally, you can find the links to the datasets below.
There are also some examples of how the CAMELS data can be used, e.g. to calculate signatures with the TOSSH toolbox (https://github.com/TOSSHtoolbox/TOSSH) and to plot them on a map.
For the plots, the BrewerMap toolbox (https://github.com/DrosteEffect/BrewerMap) is used.

To limit file size, only one P and PET product is loaded:
- for CAMELS US, PET is calculated using a calibrated Priestley-Taylor coefficient which is now set to the default value of 1.26. Observed time series are extracted from the model output files,
- for CAMELS CL, MSWEP P and HARGREAVES PET are used,
- for CAMELS GB, PET without interception correction is used,
- for CAMELS BR, CHIRPS P, GLEAM PET, and CPC T are used,
- for CAMELS AUS, AWAP P, Morton SILO PET, and AWAP T are used.

Note that some updates made are not compatible with previous versions of this repository.

Matlab 2020a was used.

## Papers describing datasets:

Addor, N., Newman, A.J., Mizukami, N. and Clark, M.P., 2017. The CAMELS data set: catchment attributes and meteorology for large-sample studies. Hydrology and Earth System Sciences (HESS), 21(10), pp.5293-5313.

Alvarez-Garreton, C., Mendoza, P.A., Boisier, J.P., Addor, N., Galleguillos, M., Zambrano-Bigiarini, M., Lara, A., Cortes, G., Garreaud, R., McPhee, J. and Ayala, A., 2018. The CAMELS-CL dataset: catchment attributes and meteorology for large sample studies-Chile dataset. Hydrology and Earth System Sciences, 22(11), pp.5817-5846.

Coxon, G., Addor, N., Bloomfield, J.P., Freer, J., Fry, M., Hannaford, J., Howden, N.J., Lane, R., Lewis, M., Robinson, E.L. and Wagener, T., 2020. CAMELS-GB: hydrometeorological time series and landscape attributes for 671 catchments in Great Britain. Earth System Science Data, 12(4), pp.2459-2483.

Chagas, V.B., Chaffe, P.L., Addor, N., Fan, F.M., Fleischmann, A.S., Paiva, R.C. and Siqueira, V.A., 2020. CAMELS-BR: hydrometeorological time series and landscape attributes for 897 catchments in Brazil. Earth System Science Data, 12(3), pp.2075-2096.

Fowler, K.J., Acharya, S.C., Addor, N., Chou, C. and Peel, M.C., 2021. CAMELS-AUS: Hydrometeorological time series and landscape attributes for 222 catchments in Australia. Earth System Science Data Discussions, pp.1-30.


## Links to datasets:

https://ral.ucar.edu/solutions/products/camels

https://doi.pangaea.de/10.1594/PANGAEA.894885?format=html#download

https://catalogue.ceh.ac.uk/documents/8344e4f3-d2ea-44f5-8afa-86d2987543a9

https://zenodo.org/record/3964745

https://doi.pangaea.de/10.1594/PANGAEA.921850
