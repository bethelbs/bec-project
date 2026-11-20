# Contributing

Thank you for your interest in this project. Questions, bug reports and suggestions are welcome as [GitHub issues](https://github.com/bethelbs/bec-project/issues).

## Reporting a problem

Please include:

- what you ran (the exact command, e.g. `pixi run reproduce`);
- what you expected and what happened instead, with the error message;
- your operating system, and the output of `pixi --version`.

## Making a change

1. Fork the repository and create a branch for your change.
2. Set up the environment with `pixi install --locked`.
3. Make the change. Keep every setting in `config/config.yaml`, every tool in `pixi.toml`, and every path relative to the project root.
4. Check that the whole analysis still runs from scratch:
   ```bash
   rm -rf data/processed results
   pixi run reproduce
   ```
5. If you change a package, commit `pixi.toml` and `pixi.lock` together.
6. If the change affects how the analysis works or what it produces, add a dated entry to `docs/decisions.md` saying what changed and why.
7. Open a pull request. The Reproduce check must pass before it is merged.

## Code style

Python scripts take their inputs, outputs and settings as command-line arguments. Run `pixi run -e pipeline snakefmt workflow/Snakefile` after editing the workflow.

Please follow the [code of conduct](CODE_OF_CONDUCT.md) in all project spaces.
