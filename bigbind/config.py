import os
from omegaconf import OmegaConf

def load_config(path):
    """ Loads configs/default.yaml, and overrides any values
    with those in the file at `path` if it exists. """

    cfg = OmegaConf.load("configs/default.yml")

    if os.path.exists(path):
        cfg = OmegaConf.merge(cfg, OmegaConf.load(path))

    # make the directories if they don't exist
    os.makedirs(cfg.data_dir, exist_ok=True)
    os.makedirs(cfg.temp_data_dir, exist_ok=True)

    # Add in some more attributes for convenience
    cfg.db_file = os.path.join(cfg.data_dir, "bigbind.db")

    return cfg

CONFIG = load_config("../configs/config.yml")