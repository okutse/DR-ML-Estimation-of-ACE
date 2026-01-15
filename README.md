This repository holds files for the Doubly-robust ML estimation of the average causal effect of differentially expressed genes (DEGs) in breast cancer


# Python Env Setup
## Install uv

```bash
# On macOS and Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# On Windows
powershell -c "irm https://astral.sh/uv/install.ps1 | iex"
```

## Setup with uv

```bash
# Sync dependencies from pyproject.toml
uv sync

# Activate the virtual environment
# On macOS/Linux
source .venv/bin/activate

# On Windows
.venv\Scripts\activate
```
