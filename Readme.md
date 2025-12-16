# Introduction
This project is used for analysis of Cross Sectional Flux method applied to TROPOMI data. 

Full report is available on the EC&V web site [Report](http:\\www.ecandv.com.au\TROPOMI-CSF)

# Project installation
Installation notes target VS Code on Windows 10 and 11. As all code is written in Python this should be easily modifed to Linux or MacOS operating system and your favourite Python IDE.

## Pre-requsites

1. Install Git. This requires administrative privilege. You will need it to download this project code.
```
winget install --id Git.Git -e --source winget
```
Restart computer.

2. From Git command configure your access to GitHub.
```
git config user.name "Jane Smith"
git config user.email janesmith@example.com
```

3. Install VS Code. This is the software development environment I use. Pick your own.
```
winget install --id Microsoft.VisualStudioCode -e --source winget
```
Restart computer.

4. Install miniconda. This is an installation software which allows compilation of some of the Python packages.
```
winget install --id Anaconda.Miniconda3 -e --source winget
```
Restart computer.

5. Install AWSCLI (AWS Command Line Interface). You will need it to download TROPOMI files from Copernicus Datastore
```
winget install --id Amazon.AWSCLI -e --source winget
```

6. Create an account with Copernicus Datastore to download ERA 5 reanalysis files [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/).

7. Using login from the previous point obtain Copernicus credentials to download TROPOMI files. Instruction on configuring access are available from [Copernicus Dataspace S3](https://documentation.dataspace.copernicus.eu/APIs/S3.html). 

## Project
1. Create a project folder for your project for example D:\Source\Repos\TROPOMI_CSF
```
mkdir D:\Source\Repos\TROPOMI_CSF
```

2. Using Visual Studio Code go to menu File -> Open Folder and open the project folder. You can to the same from a command line
```
code D:\Source\Repos\TROPOMI_CSF
```

3. In VS Code go to menu Terminal -> New Terminal. In terminal window initialize Git and download the code into the project folder.
```
git init
git remote add origin https://github.com/ECANDV/TROPOMI_CSF.git
git fetch origin master
```

4. Open VS Code Terminal and check if conda environment works.
```
conda info --envs
```

If you get an error that conda is not found manually add Conda environemnt variables to setting.json file located in .vscode folder. This assumes that conda is installed in your directory %USERPROFILE%\miniconda3

```
  "terminal.integrated.env.windows": {
        "CONDA_DEFAULT_ENV": "base",
        "CONDA_EXE": "~\\miniconda3\\Scripts\\conda.exe",
        "CONDA_PREFIX": "~\\miniconda3",
        "CONDA_PROMPT_MODIFIER": "(base)",
        "CONDA_PYTHON_EXE": "~\\miniconda3\\python.exe",
        "CONDA_SHLVL": "1",
        "PATH": "~\\miniconda3;~\\miniconda3\\Library\\mingw-w64\\bin;~\\miniconda3\\Library\\usr\\bin;~\\miniconda3\\Library\\bin;~\\miniconda3\\Scripts;~\\miniconda3\\bin;~\\miniconda3\\condabin;${env:PATH}"
  }
```
Close and re-open all VS Code Terminals as this change is only loaded when you open a new terminal.

5. In VS Code Terminal create and update environemnt, install all required projects
```
conda update conda
conda create -n TROPOMI_CSF
conda activate TROPOMI_CSF
conda install -n TROPOMI_CSF conda-forge::boto3
conda install -n TROPOMI_CSF conda-forge::gdal
conda install -n TROPOMI_CSF conda-forge::cartopy
conda install -n TROPOMI_CSF conda-forge::netCDF4
conda install -n TROPOMI_CSF conda-forge::cdsapi
```

6. In Visual Studio Code set interpreter. Enter Crtl-Shft-P, select "Python: Select Interpreter" and set to Python supplied by Miniconda "~\miniconda3\env\TROPOMI_CSF\Python.exe".  
Close and restart any VS Code Terminals.

7. Confirm that the proper environment has been set and packages installed. In VS Code terminal use the following commands. You should see * next to environemnt called TROPOMI_CSF.
```
conda info --envs
conda search --envs boto3
conda search --envs gdal
conda search --envs cartopy
conda search --envs netcdf4
conda search --envs cdsapi
```

8. Install VS Code extensions. I use the following extensions Pylance, Python, Python Debugger, Rainbow CSV, Visual Studio Keymap, but this is just how like to work. To do so from the VS Code terminal execute the following commands:

```
code --install-extension ms-python.vscode-pylance
code --install-extension ms-python.python
code --install-extension ms-python.debugpy
code --install-extension mechatroner.rainbow-csv
code --install-extension ms-vscode.vs-keybindings
```

9. Tests installation
To confirm installation you can run a sequence of tests from the VS Code terminal. The following list provides recommended order of test execution as there are dependencies between some of tests. 

Test ERA will download ERA files for case HailCreek 09956. Output of this test is used in tests: Meteorology, Algorithm_CSF. 
Output of this test is used in tests: TROPOMI, Algorithm_CSF. This can take a few minutes.

```
python -m unittest Chemistry_test.py  
python -m unittest ERA5_test.py  
python -m unittest Meteorology_test.py  
python -m unittest TROPOMI_test.py  
```

If all test pass you can display TROPOMI chart for Hail Creek orbit 09956, run CSF algorithm and list its results, using the following commands:
```
python .\TROPOMI.py chart domain HailCreek 09956
python .\Algorithm_CSF.py run HailCreek 09956 
python .\Algorithm_CSF.py list HailCreek 09956 
```

# Folder Structure
This project uses specific organization of directories. Location of these can be configured in the Config.py file
```
    Root_folder = path.join("F:")
    BoM_folder = path.join(Root_folder, "BoM")
    CSF_folder = path.join(Root_folder, "CSF")
    ERA5_folder = path.join(Root_folder, "ERA5")
    HYSPLIT_folder = path.join(Root_folder, "HYSPLIT")
    LOG_folder = path.join("Log")
    TROPOMI_folder = path.join(Root_folder, "TROPOMI")
    Paper_folder = path.join("F:", "Paper_Original")
    TROPOMI_Processors = ["010202", "010300", "010301", "010302"] 
```

Each directory contains data related to specific data source.
Log directory contains processing logs
Paper directory contains outputs used for data and charts contained in the associated report.
TROPOMI_Processors defines an array of TROPOMI processor that will be used during algorithm processing 

# Code Structure
Files can be split into 4 different categories: Programs, Metadata and Utilities, Configuration, Documentation and Tests.

## Programs
Programs can be executed from command line. Adding option -h will provide description of options for each program. For example:

```
python .\Algorithm_CSF.py -h
```

Algorithm_CSF.py        - Cross sectional flux algorithm aiming at replicating results published by Sadavarte et al 2021 \
Algorithm_HYSPLIT.py    - HYSPLIT trajectory data \
Algorithm_IME.py        - Integrated mass enhancement algorithm using CSF data \
Algorithm_TM.py         - Total mass algorithm using CSF and HYSPLIT trajectory data \
BOM_AWS.py              - Management for AWS data. You will need to buy your data from BoM as they are licensed. \
BOM_MSLP.py             - Download of publicly available MSLP data from BoM archive \
Copernicus_Catalogue.py - Search for TROPOMI files in the Copernicus Data Store \
Copernicus_S3.py        - Download TROPOMI files from the Copernicus Data Store \
ERA5.py                 - ERA data management \
Paper.py                - Chart and listing functions for charts and tables used in the report. \
SRTM.py                 - Component managing topography data coming from Shuttle Radar Topography Mission. \
TROPOMI.py              - Manage TROPOMI data \
TROPOMI_Filter.py       - Algorithm aiming at replicating filtering of orbits as published by Sadavarte et al 2021

## Metadata and Utilities
AreaCalculator.py       - Calculate area of polygons in m2 using geodetic coordinates \
Chemistry.py            - Utility converting dry column ppb to mass \
Constants.py            - Set of constants used by calculations \
Geometry.py             - Geometrical relationship between domain and pixels \
Mask.py                 - Masking for enhancement \
Meteorology.py          - Calculations of pressure averaged boundary layer wind \
Source.py               - Definition of sources of emissions. Currently contains only Hail Creek \
trajectory.py           - Utility to draw trajectories

## Configuration
Config.py               - Manage configuration of project structure  

## Documentation
LICENSE.md              - MIT license for this software \
Readme.md               - This readme file

# Licenses
This project is licensed under terms of the MIT license. I exstensively depends on Python Standard Library, Additionally it depends on the following packages which carry spearate licences.

Boto3 - licence is Apaches 2.0, available at [Boto3 GitHub](https://github.com/boto/boto3) \
Cartopy  - licence is BSD 3-Clause License, available at [Cartopy Licence](https://scitools.org.uk/cartopy/docs/latest/copyright.html) \
cdsapi - licence is Apache 2.0, available at [cdsapi GitHub](https://github.com/ecmwf/cdsapi/blob/master/LICENSE.txt) \
GDAL - licence is MIT, available at [GDAL Licence](https://gdal.org/en/stable/license.html)
Some files are licensed under BSD 2-clause, BSD 3-clause or other non-copyleft licenses. The full terms are available in the LICENSE.TXT file.
Note however that GDAL can be built against many dependencies, each of them with their own licensing terms (possibly LGPL, GPL or proprietary), hence the use of a GDAL binary can be subject to less permissive licensing terms than MIT, and it is the responsibility of users to check that they comply with the overall licensing terms. \
matplotlib - maintains its own licence whicch is available at [matplotlib Licence](https://matplotlib.org/stable/project/license.html) \
netCDF4-python - licence is MIT, available at [netCDF4-python GitHub](https://github.com/Unidata/netcdf4-python/blob/master/LICENSE) \
numpy - uses custom licence, available at [numpy GitHub](https://github.com/numpy/numpy/blob/main/LICENSE.txt) \
osgeo - licence is BSD-3, available at [osgeo Licence](https://www.osgeo.org/about/licenses/) \
PIL - licence is MIT-CMU, available at [Pillow Licence](https://raw.githubusercontent.com/python-pillow/Pillow/main/LICENSE) \
pyproj - licence is MIT, available at [PyProj GitHub](https://github.com/pyproj4/pyproj/blob/main/LICENSE) \
pystac  - licence is Apache 2.0, available at home page [pystac GitHub](https://github.com/stac-utils/pystac) \
requests - licence is Apache 2.0, available at [requests GitHub](https://github.com/psf/requests/blob/main/LICENSE)  \
shapely - licence is BSD-3, available at [shapely GitHub](https://github.com/shapely/shapely/blob/main/LICENSE.txt)