# tables

Tables! We all love them and couldn't function without them in our daily lives. Now sourmash has tables too!

## Installation

#### For regular install:

```{python}
pip install sourmash_plugin_tables
```

#### For conda environment:

```{python}
conda activate <env-name>

python -m pip install sourmash_plugin_tables
```

Or, use a conda environment yaml file:

```{yaml}
name: stable
channels:
channels:
    - conda-forge
    - bioconda
    - defaults
dependencies:
    - python>=3.10,<3.12
    - sourmash>=4.8.11,<5
    - polars
    - pip
    - pip:
        - sourmash_plugin_tables
```

#### Sanity check

```{bash}
sourmash info -v
```

Expected output:

```{bash}
== This is sourmash version 4.9.4. ==
== Please cite Irber et. al (2024), doi:10.21105/joss.06830. ==

sourmash version 4.9.4
- loaded from path: /home/colton/miniconda3/envs/smash/lib/python3.14/site-packages/sourmash/cli

khmer version: None (internal Nodegraph)

screed version 1.1.3
- loaded from path: /home/colton/miniconda3/envs/smash/lib/python3.14/site-packages/screed

the following plugins are installed:

plugin type          from python module             v     entry point name    
-------------------- ------------------------------ ----- --------------------
sourmash.cli_script  sourmash_plugin_tables         0.7   compare_rows        
sourmash.cli_script  sourmash_plugin_tables         0.7   gather_tables       
sourmash.cli_script  sourmash_plugin_tables         0.7   hash_tables         
sourmash.cli_script  sourmash_plugin_tables         0.7   prefetch_tables     
```

## Usage



## Support

We suggest filing issues in [the main sourmash issue tracker](https://github.com/dib-lab/sourmash/issues) as that receives more attention!

## Dev docs

`tables` was developed from https://github.com/sourmash-bio/sourmash_plugin_template.

### Testing

Run:
```
pytest tests
```

Example/test data located in `tests/test-data` (shocking, I know)

### Generating a release

Bump version number in `pyproject.toml` and push.

Make a new release on github.

Then pull, and:

```
make dist
```

followed by `twine upload dist/...`.
