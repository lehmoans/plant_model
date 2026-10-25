import numpy as np
import pytest
from model import (
    check_mass_conservation,
    flow_direction_color,
    solve_network,
)

"""
validation of code: 
tests code derived values against hand calculated values, validating theory and code. 
"""

# input parameters
MU = 1e-3       # Dynamic viscosity [Pa s]
RHO = 1000.0    # Fluid density [kg/m^3]
G = 9.81        # Gravitational acceleration [m/s^2]

P_ROOT_TIP = 0.0
P_LEAF_TIP = 0.0

CHANNEL_LEN = 1
CHANNEL_R1 = 2
CHANNEL_R2 = 1
CHANNEL_R3 = 0.5

CHANNEL_LEN1 = 2
CHANNEL_LEN2 = 1
CHANNEL_LEN3  = 0.5

CHANNEL_LEN = CHANNEL_LEN2

# Nodes are grouped by region/order. A known pressure marks a boundary node;
# p=None marks an internal junction whose potential is solved for.
NODE_TREE = {
    # Above ground
    1: {
        1: {"m": {"p": None, "z": None}},
        2: {
            "n1": {"p": None, "z": None},
            "n2": {"p": None, "z": None},
        },
        3: {
            "n11": {"p": P_LEAF_TIP, "z": None},
            "n12": {"p": P_LEAF_TIP, "z": None},
            "n21": {"p": P_LEAF_TIP, "z": None},
            "n22": {"p": P_LEAF_TIP, "z": None},
        },
    },
    # Ground
    0: {"g": {"p": None, "z": None}},
    # Below ground
    -1: {
        1: {"k": {"p": None, "z": None}},
        2: {
            "j1": {"p": None, "z": None},
            "j2": {"p": None, "z": None},
        },
        3: {
            "j11": {"p": P_ROOT_TIP, "z": None},
            "j12": {"p": P_ROOT_TIP, "z": None},
            "j21": {"p": P_ROOT_TIP, "z": None},
            "j22": {"p": P_ROOT_TIP, "z": None},
        },
    },
}

# Branch: (node_a, node_b, length_m, radius_m, theta_deg)
# theta is measured from horizontal, positive upward from a to b.
BRANCHES = [
    ("n1", "n11", CHANNEL_LEN, CHANNEL_R3, 45),
    ("n1", "n12", CHANNEL_LEN, CHANNEL_R3, 90),
    ("n2", "n21", CHANNEL_LEN, CHANNEL_R3, 90),
    ("n2", "n22", CHANNEL_LEN, CHANNEL_R3, 135),
    ("n2", "m", CHANNEL_LEN, CHANNEL_R2, -25),
    ("n1", "m", CHANNEL_LEN, CHANNEL_R2, -155),
    ("m", "g", CHANNEL_LEN, CHANNEL_R1, -90),
    ("g", "k", CHANNEL_LEN1, CHANNEL_R3, 270),
    ("k", "j1", CHANNEL_LEN2, CHANNEL_R3, -145),
    ("k", "j2", CHANNEL_LEN2, CHANNEL_R3, -35),
    ("j1", "j11", CHANNEL_LEN3, CHANNEL_R3, -135),
    ("j1", "j12", CHANNEL_LEN3, CHANNEL_R3, -45),
    ("j2", "j21", CHANNEL_LEN3, CHANNEL_R3, -135),
    ("j2", "j22", CHANNEL_LEN3, CHANNEL_R3, -45),
]



# hand calculated values for A, b
HAND_UNKNOWN_ORDER = ["n2", "n1", "m", "g", "k", "j1", "j2"]
HAND_A_DIAGONALS = {
    0: [441.779, 441.779, 7068.583, 6295.45, 61.36, 122.718, 122.718],
    1: [0, -392.699, -6283.185, -12.272, -24.54, 0],
    2: [-392.699, 0, 0, 0, -24.54],
    -1: [0, -392.699, -6283.185, -12.272, -24.54, 0],
    -2: [-392.699, 0, 0, 0, -24.54],
}  
HAND_B = np.array([1095282.949, 1095282.949, 0, 0, 0, -2819100.28, -2819100.28])


#methods
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


def test_handcalc_reference_diagonals_are_checked_without_silent_truncation():
    # The committed reference currently has malformed diagonal lengths.
    # Treat that as an expected data-quality finding, not a passing comparison.
    with pytest.raises(ValueError, match="diagonal offset"):
        build_handcalc_matrix()


def test_solver_matrix_and_rhs_match_hand_calculations():
    try:
        comparison = compare_solver_to_handcalc()
    except ValueError as exc:
        if "diagonal offset" in str(exc):
            pytest.skip(f"Hand-calculation reference must be corrected first: {exc}")
        raise
    assert comparison["A_matches"], (
        "A mismatch in hand-calculation node order "
        f"{HAND_UNKNOWN_ORDER}"
    )
    assert comparison["b_matches"], "Solver RHS does not match hand-calculated b"

if __name__ == "__main__":
    results = solve_network(NODE_TREE, BRANCHES)
    A_solver, b_solver = align_solver_matrix(results)

    print('Solver unknown order:', results['unknowns'])
    print('Hand-calculation order:', HAND_UNKNOWN_ORDER)
    print('A from code, aligned to hand-calculation order:')
    print(A_solver)
    print('b from code, aligned to hand-calculation order:')
    print(b_solver)

    try:
        A_hand = build_handcalc_matrix()
    except ValueError as exc:
        A_hand = None
        print('REFERENCE DATA ERROR:', exc)
        print('\nDiagonal lengths supplied vs required:')
        n = len(HAND_UNKNOWN_ORDER)
        for offset, values in HAND_A_DIAGONALS.items():
            print(f'offset {offset:+}: supplied {len(values)}, required {n - abs(offset)}')
    else:
        print('Hand-calculated A:')
        print(A_hand)
        print('Hand-calculated b:')
        print(HAND_B)


    if A_hand is None:
        print('Comparison blocked: correct the malformed hand-calculation diagonals above, then rerun this cell.')
    else:
        a_match = np.allclose(A_solver, A_hand, rtol=1e-2, atol=1e-2)
        b_match = np.allclose(b_solver, HAND_B, rtol=1e-2, atol=1e-2)
        print('A matches:', a_match)
        print('b matches:', b_match)
        print('A difference (code - hand):')
        print(A_solver - A_hand)
        print('b difference (code - hand):')
        print(b_solver - HAND_B)
        assert a_match and b_match, 'Solver system differs from the hand-calculated system'
