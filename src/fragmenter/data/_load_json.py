import json
import os
from pathlib import Path
from typing import Union

from fragmenter.typing import SMARTS

FILEPATH_JSON_DATA = Path(__file__).parent.resolve() / "json"
if not FILEPATH_JSON_DATA.is_dir():
    raise ValueError("'data' is not packaged")


def __load_json_atomic_groups(
    filepath_json: Union[os.PathLike, str],
) -> dict[str, list[SMARTS]]:
    with open(filepath_json, mode="r", encoding="utf-8") as f:
        return json.load(f)


SMARTS_UNIFAC = __load_json_atomic_groups(
    FILEPATH_JSON_DATA / "UNIFAC-ATOMIC-GROUPS.json"
)
"""UNIFAC's atomic groups

- key: name of atomic group
- value: SMARTS patterns

"""
