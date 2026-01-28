"""Shared utilities for command modules."""


def parse_value(value_str):
    """Parse a string value into appropriate Python type."""
    if value_str.lower() == "true":
        return True
    if value_str.lower() == "false":
        return False
    if value_str.lower() == "none":
        return None

    try:
        return int(value_str)
    except ValueError:
        pass

    try:
        return float(value_str)
    except ValueError:
        pass

    return value_str


def parse_extra_args(extra_args):
    """Parse a list of extra arguments into kwargs dict."""
    kwargs = {}
    i = 0
    while i < len(extra_args):
        arg = extra_args[i]
        if arg.startswith("--"):
            key = arg[2:].replace("-", "_")
            if i + 1 < len(extra_args) and not extra_args[i + 1].startswith("--"):
                kwargs[key] = parse_value(extra_args[i + 1])
                i += 2
            else:
                kwargs[key] = True
                i += 1
        else:
            i += 1
    return kwargs
