==========
gridmetETL
==========


.. image:: https://img.shields.io/pypi/v/gridmetetl.svg
        :target: https://pypi.python.org/pypi/gridmetetl

.. image:: https://img.shields.io/travis/rmcd-mscb/gridmetetl.svg
        :target: https://travis-ci.com/rmcd-mscb/gridmetetl

.. image:: https://readthedocs.org/projects/gridmetetl/badge/?version=latest
        :target: https://gridmetetl.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




Extract gridmet, Translate to hru, Load to netCDF


* Free software: MIT license
* Documentation: https://gridmetetl.readthedocs.io.

Install
-------
1. Create conda env as follows
    * conda create -n gmetl python=3.7
    * conda activate gmetl
    * conda install -c conda-forge numpy matplotlib pandas geopandas xarray netcdf4 requests dask
    * conda install -c conda-forge jupyterlab
   OPTIONAL:
    * conda install -c conda-forge git
    * conda install -c conda-forge pip

2. Clone repository
    * cd gridmetetl
    * Develop code:
        * pip install -e .
    * Use code
        * pip install .
Use
-------
From your conda environment created above:

```
(gmetl) B:\gitbmi\gridmetetl>gridmetetl -h

usage: gridmet_etl [-h] -t extraction type [-p YYYY-MM-DD) (YYYY-MM-DD]
                   [-d numdays] [-f output_file_prefix] -i input_path -o
                   output_path -w weight_file
                   [-v [GridMet_Variables [GridMet_Variables ...]]]

map gridded climate data to polygon using zonal area weighted mean

optional arguments:
  -h, --help            show this help message and exit
  -t extraction type, --extract_type extraction type
                        extract method: (days) or (date)
  -p (YYYY-MM-DD) (YYYY-MM-DD), --period (YYYY-MM-DD) (YYYY-MM-DD)
                        option: start date and end date of retrieval (YYYY-MM-
                        DD)
  -d numdays, --days numdays
                        option: number of days to retrieve; if specified take
                        precedence over -s & -e option
  -f output_file_prefix, --file_prefix output_file_prefix
                        option: prefix for output files
  -i input_path, --inpath input_path
                        input_path (location of HRU shapefiles)
  -o output_path, --outpath output_path
                        Output path (location of netcdf output files by
                        shapefile output)
  -w weight_file, --weightsfile weight_file
                        path/weight.csv - path/name of weight file
  -v [GridMet_Variables [GridMet_Variables ...]], --variables [GridMet_Variables [GridMet_Variables ...]]
                        over-ride default vars
 ```
### Do an ETL:

```
gridmetetl -t date -p 2018-09-01 2018-09-02 -i ../../GitRepos/onhm-fetcher-parser/Data -o ../../GitRepos/onhm-fetcher-parser/Output -w ../../onhm-fetcher-parser/Data/weights.csv
```

### Additional examples:
https://github.com/nhm-usgs/gridmetetl/blob/master/Examples/Example_code_usage.ipynb

Features
--------

* TODO

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
