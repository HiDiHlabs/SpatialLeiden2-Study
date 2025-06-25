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


# Plotting utils
mm = 0.1 / 2.54

fontsize = 6
rc_params = {
    "axes.titlesize": fontsize,
    "axes.labelsize": fontsize,
    "xtick.labelsize": fontsize,
    "ytick.labelsize": fontsize,
    "legend.fontsize": fontsize,
    "legend.title_fontsize": fontsize,
}
