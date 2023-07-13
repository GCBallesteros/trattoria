# How to develop Trattoria

Trattoria involves 3 repos and some compiled code so getting up and running a dev
environment for it is not as trivial as I would like it.

The following instructions are a first step towards making it easier for newcomers
and future me to setup.

1. Clone `trattoria`, `trattoria-core` and `tttr-toolbox`.
2. In `trattoria-core` modify the Cargo.toml file so that it picks up `tttr-toolbox`
   not from crates.io but from a local folder on your computer. This will have to be
   undone before you push a new version to crates.io.
3. Create a Python venv and move into it.
4. Assuming you are still inside the `trattoria-core` you should now run `poetry install`.
   This will install the required Python dependencies.
5. Now you can run `poetry run maturin develop` which will compile the rust components
   and make them available from within the venv. You can test this by opening a Python
   interpreter and `import trattoria_core`.
6. Without leaving the venv move to the `trattoria` repo. Open the `pyproject.toml`
   file and comment out the `trattoria-core` dependency.  This will temporatily disable
   the dependency on trattoria-core so that the next step doesnt attempt to fetch `trattoria-core`
   from pypi which we haven't up.
7. Now `poetry install` and finally check that everythin works by starting a REPL and
   and verifying you can `import trattoria`.

Reached this point you can start developing `trattoria-core` and `tttr-toolbox` and
test the code by calling it from Python using `trattoria`. Remember that between
changes you will have to rerun step 5 above.
