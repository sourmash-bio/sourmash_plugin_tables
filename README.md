# sourmash_plugin_tables

## Installation

```
pip install sourmash_plugin_tables
```

## Usage

non-xyz info goes here!

## Support

We suggest filing issues in [the main sourmash issue tracker](https://github.com/dib-lab/sourmash/issues) as that receives more attention!

## Dev docs

`tables` was developed from https://github.com/sourmash-bio/sourmash_plugin_template.

### Testing

Run:
```
pytest tests
```

### Generating a release

Bump version number in `pyproject.toml` and push.

Make a new release on github.

Then pull, and:

```
make dist
```

followed by `twine upload dist/...`.
