# PlioMIP2-AMOC
Codes used for analysis of the PlioMIP2 AMOC in Weiffenbach et al. (2022, in prep.)

This code has been used to process PlioMIP2 data and create the figures of the above manuscript, including the supplementary material. 

Data provided by individual PlioMIP2 modelling groups includes the ocean potential temperature, ocean meridional velocity, ocean salinity, the Atlantic meridional streamfunction, the total Atlantic OHT and the atmospheric zonal and meridional wind at 1000 hPa.
This data is available upon request from Alan M. Haywood (a.m.haywood@leeds.ac.uk), with the exception of IPSL-CM6A, EC-Earth3-LR and GISS2.1G. PlioMIP2 data from IPSL-CM6A, EC-Earth3-LR and GISS2.1G can be obtained from the Earth System Grid Federation (ESGF) (https://esgf-node.llnl.gov/search/cmip6/).
The atmospheric surface freshwater flux (precipitation minus evaporation; PmE) fields for all models except MIROC4m have been provided by Han et al. (2021).
The links to SST reconstruction data are provided in the code and manuscript. 

The "Processing" folder contains all processing that has been done on raw downloaded data. It contains three subfolders: pre-processing, interpolate_3D and freshwater_transport. 
The notebooks that are in the "Processing" folder itself contain the main processing such as the OHT separation, interpolation of surface fields and computing Atlantic zonal mean fields.
Pre-processing is done for some models to compress multiple files containing data into a single file (often quite inefficiently, but it works). 
The notebooks in the "interpolate_3D" subfolder compute 3-D salinity, velocity and potential temperature fields that are horizontally interpolated to a regular 1x1 degree grid. Only the interpolated 3-D velocity field is used in analysis (FigureS06.ipynb).
The notebooks in "freshwater_transport" subfolder are used calculate the freshwater transport components in Atlantic Ocean, as well as the freshwater transports into the Labrador Sea, across the Greenland-Scotland Ridge and through the Bering Strait.

The notebooks in the "Figures" folder contain all the figures and data that are presented in the manuscript. The filename lists the figure(s) made with the notebook.
