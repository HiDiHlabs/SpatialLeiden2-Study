from pathlib import Path

SHELL = "/bin/bash"
CONDA_PATH = "/opt/conda/bin/activate"

DATA_DIR = Path("data")
RESULT_DIR = Path("results")


def activate_env_cmd(env: str) -> str:
    return f"source {CONDA_PATH} {env}"


def remove_tick_and_label(ax):
    ax.set(xticklabels=[], yticklabels=[], xlabel=None, ylabel=None)
    ax.tick_params(left=False, bottom=False)
