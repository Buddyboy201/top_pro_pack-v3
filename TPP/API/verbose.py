
global VERBOSE
VERBOSE = False

def toggle_io():
    globals()["VERBOSE"] = not globals()["VERBOSE"]

def handle_debug(f, *args, **kwargs):
    if VERBOSE:
        f(*args, **kwargs)