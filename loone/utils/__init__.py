import os
from yaml import safe_load


def load_config(workspace: str) -> dict:
    """Load the configuration file.

    Args:
        workspace (str): Path to directory containing the config file and its associated input data.

    Returns:
        dict: dictionary containing the configuration parameters.

    Raises:
        FileNotFoundError: If the configuration file does not exist.
    """
    for config_file_name in ["config.yaml", "config.yml"]:
        config_file = os.path.join(workspace, config_file_name)
        if os.path.exists(config_file):
            with open(config_file, "r") as f:
                config = safe_load(f)
            return config
    else:
        raise FileNotFoundError("Config file not found in the workspace.")
