import os
import json


def load_symbols():
    """Load the element symbols from the json file."""
    config_dir = os.path.join(os.path.dirname(__file__), "../config")
    json_path = os.path.join(config_dir, "elements.json")

    with open(json_path, "r") as file:
        element_symbols = json.load(file)

    return element_symbols
