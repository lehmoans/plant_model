import numpy as np

from model import (
    BRANCHES,
    NODE_TREE,
    check_mass_conservation,
    flow_direction_color,
    solve_network,
)


def test_nested_node_structure_solves_and_conserves_mass():
    results = solve_network(NODE_TREE, BRANCHES)
    report = check_mass_conservation(
        results["nodes"], results["branches"], results["flows"]
    )

    internal_nodes = [item for item in report if item["is_internal"]]
    assert internal_nodes
    assert all(item["passed"] for item in internal_nodes)
    assert np.all(np.isfinite(results["flows"]))


def test_flow_plot_colours_actual_vertical_direction():
    positions = {
        "lower": (0.0, 0.0),
        "upper": (0.0, 1.0),
        "level_a": (0.0, 0.5),
        "level_b": (1.0, 0.5),
    }

    assert flow_direction_color("lower", "upper", 1.0, positions) == "tab:green"
    assert flow_direction_color("lower", "upper", -1.0, positions) == "tab:red"
    assert flow_direction_color("level_a", "level_b", 1.0, positions) == "0.55"
