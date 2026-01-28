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
    """Parse a list of extra arguments into kwargs dict.

    Supports list-like inputs: --key val1 val2 val3 becomes key=[val1, val2, val3].
    Single values are passed as-is (not wrapped in a list).
    """
    kwargs = {}
    i = 0
    while i < len(extra_args):
        arg = extra_args[i]
        if arg.startswith("--"):
            key = arg[2:].replace("-", "_")
            i += 1
            # Collect all values until the next --flag or end of args
            values = []
            while i < len(extra_args) and not extra_args[i].startswith("--"):
                values.append(parse_value(extra_args[i]))
                i += 1
            # Assign: single value as-is, multiple as list, none as True (flag)
            if len(values) == 0:
                kwargs[key] = True
            elif len(values) == 1:
                kwargs[key] = values[0]
            else:
                kwargs[key] = values
        else:
            i += 1
    return kwargs
