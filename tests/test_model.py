import numpy as np
import pytest

from model import (
    NODE_TREE,
    BRANCHES,
    check_mass_conservation,
    flow_direction_color,
    solve_network,
)

# Transcribed from the hand-calculation data committed on 2026-10-22.
# The node order is deliberately explicit so comparisons do not depend on
# the dictionary traversal order used by the solver.
HAND_UNKNOWN_ORDER = ["n2", "n1", "m", "g", "k", "j1", "j2"]
HAND_A_DIAGONALS = {
    0: [0.084, 0.084, 0.0817, 0.0811, 0.163, 0.081, 0.081],
    1: [0, -0.0025, 0, -0.081, 0, -0.041, 0],
    2: [-0.0025, 0, 0, 0, -0.041],
    -1: [0, -0.0025, 0, -0.000159, 0, -0.041, 0],
    -2: [-0.0025, 0, 0, 0, -0.041],
}
HAND_B = np.array([1785.302, 1785.302, 0, 0, 0, -1148.614, -1148.614])


def build_handcalc_matrix():
    """Build A only when each supplied diagonal has the correct length.

    Raises ValueError rather than silently truncating or padding the recorded
    hand-calculation values.
    """
    n = len(HAND_UNKNOWN_ORDER)
    A = np.zeros((n, n), dtype=float)
    for offset, values in HAND_A_DIAGONALS.items():
        expected_length = n - abs(offset)
        if len(values) != expected_length:
            raise ValueError(
                f"Hand-calculation diagonal offset {offset} has "
                f"{len(values)} values; a {n}x{n} matrix requires "
                f"{expected_length}. Please correct the transcribed reference data."
            )
        A[np.diag_indices_from(A)] = A[np.diag_indices_from(A)]  # explicit no-op
        A[np.eye(n, k=offset, dtype=bool)] = values
    return A


def align_solver_matrix(results):
    """Reorder the solver matrix/RHS to HAND_UNKNOWN_ORDER."""
    code_order = results["unknowns"]
    missing = set(HAND_UNKNOWN_ORDER) - set(code_order)
    extra = set(code_order) - set(HAND_UNKNOWN_ORDER)
    if missing or extra:
        raise AssertionError(f"Unknown-node mismatch; missing={missing}, extra={extra}")
    permutation = [code_order.index(node) for node in HAND_UNKNOWN_ORDER]
    return (
        results["matrix"][np.ix_(permutation, permutation)],
        results["rhs"][permutation],
    )


def compare_solver_to_handcalc():
    """Return the aligned solver system and hand-calculation comparison."""
    results = solve_network(NODE_TREE, BRANCHES)
    A_solver, b_solver = align_solver_matrix(results)
    A_hand = build_handcalc_matrix()
    return {
        "results": results,
        "A_solver": A_solver,
        "b_solver": b_solver,
        "A_hand": A_hand,
        "b_hand": HAND_B.copy(),
        "A_matches": np.allclose(A_solver, A_hand, rtol=1e-6, atol=1e-9),
        "b_matches": np.allclose(b_solver, HAND_B, rtol=1e-6, atol=1e-6),
    }


def test_nested_node_structure_solves_and_conserves_mass():
    results = solve_network(NODE_TREE, BRANCHES)
    report = check_mass_conservation(results["nodes"], results["branches"], results["flows"])
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


def test_handcalc_reference_diagonals_have_valid_lengths():
    # This intentionally fails while the committed reference data is malformed.
    # It prevents a misleading comparison made by silently dropping values.
    build_handcalc_matrix()


def test_solver_matrix_and_rhs_match_hand_calculations():
    comparison = compare_solver_to_handcalc()
    assert comparison["A_matches"], (
        "A mismatch in hand-calculation node order "
        f"{HAND_UNKNOWN_ORDER}"
    )
    assert comparison["b_matches"], "Solver RHS does not match hand-calculated b"
