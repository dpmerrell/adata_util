"""Tests for shared utility functions."""

from adata_util.commands.utils import parse_value, parse_extra_args


def test_parse_value_bool():
    """Test parsing boolean values."""
    assert parse_value("true") is True
    assert parse_value("True") is True
    assert parse_value("false") is False
    assert parse_value("False") is False


def test_parse_value_none():
    """Test parsing None values."""
    assert parse_value("none") is None
    assert parse_value("None") is None


def test_parse_value_int():
    """Test parsing integer values."""
    assert parse_value("42") == 42
    assert parse_value("-10") == -10
    assert parse_value("0") == 0


def test_parse_value_float():
    """Test parsing float values."""
    assert parse_value("3.14") == 3.14
    assert parse_value("-0.5") == -0.5


def test_parse_value_string():
    """Test parsing string values."""
    assert parse_value("hello") == "hello"
    assert parse_value("some_value") == "some_value"


def test_parse_extra_args_single_value():
    """Test parsing a single value for a key."""
    result = parse_extra_args(["--color", "cluster"])
    assert result == {"color": "cluster"}


def test_parse_extra_args_multiple_values():
    """Test parsing multiple values as a list."""
    result = parse_extra_args(["--color", "cluster", "cell_type"])
    assert result == {"color": ["cluster", "cell_type"]}


def test_parse_extra_args_flag():
    """Test parsing a flag with no value."""
    result = parse_extra_args(["--show"])
    assert result == {"show": True}


def test_parse_extra_args_mixed():
    """Test parsing mixed arguments."""
    result = parse_extra_args([
        "--color", "cluster", "cell_type",
        "--n_comps", "10",
        "--show",
    ])
    assert result == {
        "color": ["cluster", "cell_type"],
        "n_comps": 10,
        "show": True,
    }


def test_parse_extra_args_hyphen_to_underscore():
    """Test that hyphens in keys are converted to underscores."""
    result = parse_extra_args(["--my-key", "value"])
    assert result == {"my_key": "value"}


def test_parse_extra_args_empty():
    """Test parsing empty args."""
    result = parse_extra_args([])
    assert result == {}
