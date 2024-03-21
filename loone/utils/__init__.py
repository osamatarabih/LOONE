from yaml import safe_load


def load_config(config_path: str) -> dict:
    """Load the configuration file.

    Args:
        config_path (str): Path to the configuration file.

    Returns:
        dict: dictionary containing the configuration parameters.

    Raises:
        FileNotFoundError: If the configuration file does not exist.
    """
    try:
        with open(config_path, "r") as f:
            config = safe_load(f)

        return config
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Error loading configuration file: {e}")
