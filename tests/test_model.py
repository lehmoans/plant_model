import numpy as np


from model import (
    check_mass_conservation,
    flow_direction_color,
    solve_network,
)
#################
#data to test code calcs
MU = 1e-3       # Dynamic viscosity [Pa s]
RHO = 1000.0    # Fluid density [kg/m^3]
G = 9.81        # Gravitational acceleration [m/s^2]

P_ROOT_TIP = 0.0
P_LEAF_TIP = 0.0

CHANNEL_LEN1 = 2
CHANNEL_LEN2 = 1
CHANNEL_LEN3  = 0.5

CHANNEL_R1 = 2
CHANNEL_R2 = 1
CHANNEL_R3 = 0.5

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
    ("g", "k", CHANNEL_LEN, CHANNEL_R, 270),
    ("k", "j1", CHANNEL_LEN, CHANNEL_R, -145),
    ("k", "j2", CHANNEL_LEN, CHANNEL_R, -35),
    ("j1", "j11", CHANNEL_LEN, CHANNEL_R, -135),
    ("j1", "j12", CHANNEL_LEN, CHANNEL_R, -45),
    ("j2", "j21", CHANNEL_LEN, CHANNEL_R, -135),
    ("j2", "j22", CHANNEL_LEN, CHANNEL_R, -45),
]

# A matrix diag values, 0 elsewhere
A_diag  = {
    0: [0.084, 0.084, 0.0817, 0.0811,0.163,0.081,0.081],
    1: [0,-0.0025,0,-0.081,0,-0.041,0],
    2:[-0.0025,0,0,0,-0.041],
    -1:[0,-0.0025,0,-0.000159,0,-0.041,0],
    -2: [-0.0025,0,0,0,-0.041],
}

b_calc = [1785.302, 1785.302,0,0,0,-1148.614,-1148.614]


################# end of data

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

def checkA_against_handcalc(A, unknowns):

    n2  = n1 = m =g = k1 = j1 = j2 = None
    phi  = [n2, n1,m,g,k1,j1,j2] #unknowns

    # A matrix diag values, 0 elsewhere
    A_diag  = {
        0: [0.084, 0.084, 0.0817, 0.0811,0.163,0.081,0.081],
        1: [0,-0.0025,0,-0.081,0,-0.041,0],
        2:[-0.0025,0,0,0,-0.041],
        -1:[0,-0.0025,0,-0.000159,0,-0.041,0],
        -2: [-0.0025,0,0,0,-0.041],
    }

    b_calc = [1785.302, 1785.302,0,0,0,-1148.614,-1148.614]
    len_n = len(phi)
    A_handcalc = np.zeros(len_n, len_n)


    for diag_offset, value in A_diag.items():
        A_handcalc[np.eye(len_n, k = diag_offset,dtype=bool)] = value

    #permutation
    code_order = unknowns
    check_order = phi

    # Indices that convert check ordering to code ordering
    perm = [check_order.index(node) for node in code_order]

    A_check_aligned = A_handcalc[np.ix_(perm, perm)]

    assert np.allclose(A, A_check_aligned), "Code does not produce hand calculated A values"



