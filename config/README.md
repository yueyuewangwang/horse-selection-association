# Configuration

```bash
cp config/config.example.env config/config.env
# edit config/config.env
source config/config.env
bash workflow/run_stage.sh --help
```

`PROJECT_ROOT` is only a namespace for the user's local data tree. Inputs may
be stored elsewhere by changing the individual variables. The private
`config/config.env` file is ignored and should not be committed.
