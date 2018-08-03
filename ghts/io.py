import os.path
import json
import numpy as np

class NumpyEncoder(json.JSONEncoder):
    """Thanks to https://stackoverflow.com/a/47626762/7909981"""

    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


class NumpyDecoder(json.JSONDecoder):
    """Thanks to https://stackoverflow.com/questions/10885238"""

    def __init__(self, list_type=list,  **kwargs):
        json.JSONDecoder.__init__(self, **kwargs)
        # Use the custom JSONArray
        self.parse_array = self.JSONArray
        self.scan_once = json.scanner.py_make_scanner(self)

    def JSONArray(self, s_and_end, scan_once, **kwargs):
        values, end = json.decoder.JSONArray(s_and_end, scan_once, **kwargs)
        return np.array(values), end


def load_state(filename):
    if not os.path.isfile(filename):
        return None
    with open(filename, 'r') as f:
        return json.load(f, cls=NumpyDecoder)


def dump_state(state, filename):
    with open(filename, 'w') as f:
        json.dump(state, f, cls=NumpyEncoder, indent=4)