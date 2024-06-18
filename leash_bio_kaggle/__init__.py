"""Awesome `leash_bio_kaggle` is a Python package created with https://gitlab.com/oasci/software/snakeframe"""

from typing import Any

import os
import sys
from ast import literal_eval

from loguru import logger

logger.disable("leash_bio_kaggle")

LOG_FORMAT = (
    "<green>{time:HH:mm:ss}</green> | "
    "<level>{level: <8}</level> | "
    "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>"
)


def enable_logging(
    level_set: int, stdout_set: bool = True, file_path: str | None = None
) -> None:
    r"""Enable logging.

    Args:
        level: Requested log level: `10` is debug, `20` is info.
        file_path: Also write logs to files here.
    """
    config: dict[str, Any] = {"handlers": []}
    if stdout_set:
        config["handlers"].append(
            {"sink": sys.stdout, "level": level_set, "format": LOG_FORMAT}
        )
    if isinstance(file_path, str):
        config["handlers"].append(
            {"sink": file_path, "level": level_set, "format": LOG_FORMAT}
        )
    # https://loguru.readthedocs.io/en/stable/api/logger.html#loguru._logger.Logger.configure
    logger.configure(**config)

    logger.enable("leash_bio_kaggle")


if literal_eval(os.environ.get("leash_bio_kaggle_LOG", "False")):
    level = int(os.environ.get("leash_bio_kaggle_LOG_LEVEL", 20))
    stdout = literal_eval(os.environ.get("leash_bio_kaggle_STDOUT", "True"))
    log_file_path = os.environ.get("leash_bio_kaggle_LOG_FILE_PATH", None)
    enable_logging(level, stdout, log_file_path)
