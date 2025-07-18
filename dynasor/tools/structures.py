import numpy as np
from ase import Atoms
from dynasor.modes.atoms import Prim


def align_structure(atoms: Atoms, atol: float = 1e-5):
    """

    Rotate and realign atoms object such that
    * the first cell vector points along the x-directon
    * the second cell vector lies in the xy-plane

    Modifies input ``atoms`` object in place.

    Parameters
    ----------
    atoms
        input structure to be rotated aligned wtih the x,y,z coordinte system
    atol
        absolute tolerance used for sanity checking the cell
    """
    _align_a_onto_xy(atoms, atol)
    _align_a_onto_x(atoms, atol)
    _align_b_onto_xy(atoms, atol)


def get_offset_index(primitive, supercell, tol=0.01, wrap=True):
    """ Get the basis index and primitive cell offsets for a supercell


    Simple iteration which should be fairly quick. If more stability is needed consider:

    P = ideal.cell @ prim.cell_inv.T find the P-matrix
    P *= len(ideal)/len(prim)/det(P) compensate for strain
    ref_atoms = make_supercell(round(P), prim)
    find assignment using ref_atoms and ideal via hungarian algorithm and mic distances


    Parameters
    ----------
    primitive
        primitive cell as either ASEAtoms
    supercell
        some ideal repetition of the prim cell and possible wrapping  as either ASEAtoms
    tol
        tolerance length parameter. Increase to allow for slgihtly rattled or strained cells.
    wrap
        It might happen that the ideal cell boundary cuts through a unit cell
        whos lattice points lie within the ideal cell. If there is a basis one
        atom belonging to this unit cell might get wrapped while the other is
        not. Then the wrapped atom now belong to a lattice point outside the P
        matrix so to say. This would result in more lattice points than
        expected from N_unit = len(ideal)/len(prim).

    Returns
    -------
    offsets
        The lattice points as integers in (N, 3)-array
    index
        The basis indices as integers in (N,)-array
    """

    if not isinstance(primitive, Atoms):
        raise ValueError
    if not isinstance(supercell, Atoms):
        raise ValueError

    prim = Prim(primitive)

    from dynasor.tools.structures import get_P_matrix
    from dynasor.modes.tools import inv

    P = get_P_matrix(primitive.cell, supercell.cell)  # P C = S
    P_inv = inv(P)

    lattice, basis = [], []
    # Pick an atom in the supercell
    for pos_ideal in supercell.positions:
        # Does this atom perhaps belong to site "index"?
        for index, pos_prim in enumerate(primitive.positions):
            # if so we can remove the basis position vector and should end up on a lattice site
            diff_pos = pos_ideal - pos_prim
            # The lattice site has integer coordinates in reduced coordinates
            prim_spos = diff_pos @ prim.inv_cell.T
            # Rounding should not affect the coordinate much if it is integer
            prim_spos_round = np.round(prim_spos).astype(int)
            # If the rounded spos and unrounded spos are the same
            if np.allclose(prim_spos, prim_spos_round, rtol=0, atol=tol):
                # Since P_inv is represented using fractions we can neatly
                # write the supercell spos of the lattice point using fractions
                # and easily determine if it needs wrapping or not without
                # worry about numerics
                ideal_spos = prim_spos_round @ P_inv
                # wrap if needed
                ideal_spos_wrap = ideal_spos % 1 if wrap else ideal_spos
                # This should be integer again
                prim_spos_wrap = (ideal_spos_wrap @ P).astype(int)
                # add the results and break out from the basis site loop
                lattice.append(prim_spos_wrap)
                basis.append(index)
                break
        else:  # we get here by not breaking out from the basis site loop.
            # This means that the candidate lattice site where not close to integers

            raise ValueError(f' {prim_spos} {prim_spos_round} Supercell not compatible '
                             'with primitive cell.')

    lattice = np.array(lattice)
    basis = np.array(basis)

    # We should have found len(ideal) unique positions
    lattice_basis = [tuple((*lp, i)) for lp, i in zip(lattice, basis)]
    assert len(set(lattice_basis)) == len(supercell)

    return lattice, basis


def get_P_matrix(c, S):
    """ P c = S  ->  c.T P.T = S.T

    The P matrix must be an integer matrix. c is the primitive cell metric and
    S is the supercell metric as row-vectors
    """
    PT = np.linalg.solve(c.T, S.T)
    P_float = PT.T
    P = np.round(P_float).astype(int)
    if not np.allclose(P_float, P) or not np.allclose(P @ c, S):
        raise ValueError(
            f'Please check that the supercell metric ({S}) is related to the'
            f' the primitive cell {c} by an integer transformation matrix.')
    return P


def _align_a_onto_xy(atoms, atol):
    """Rotate cell so that a is in the xy-plane"""

    # get angle towards xy
    # will break if a is along z
    assert np.any(atoms.cell[0, :2])

    cell = atoms.cell.array.copy()

    a = cell[0]
    a_xy = a.copy()
    a_xy[2] = 0  # projection of a onto xy-plane

    # angle between a and xy-plane
    cosa = np.dot(a, a_xy) / np.linalg.norm(a) / np.linalg.norm(a_xy)

    # cosa should be in the interval (0, 1]
    assert not np.isclose(cosa, 0)
    if cosa > 1:
        assert np.isclose(cosa, 1)
    cosa = min(cosa, 1)
    cosa = max(cosa, 0)

    # angle between a and xy-plane in degs
    angle_xy_deg = np.rad2deg(np.arccos(cosa))

    # get unit vector to rotate around
    vec = np.cross(a_xy, [0, 0, 1])
    vec = vec / np.linalg.norm(vec)
    assert vec[2] == 0

    # Determine if the rotation should be positive or negative depending on
    # whether a is pointing in the +z or -z direction
    sign = -1 if a[2] > 0 else +1

    # rotate
    atoms.rotate(sign * angle_xy_deg, vec, rotate_cell=True)

    assert np.isclose(atoms.cell[0, 2], 0, atol=atol, rtol=0), atoms.cell


def _align_a_onto_x(atoms, atol):
    assert np.isclose(atoms.cell[0, 2], 0, atol=atol, rtol=0)  # make sure a is in xy-plane

    a = atoms.cell[0]
    a_x = a[0]
    a_y = a[1]

    # angle between a and x-axis (a is already in xy-plane)

    # tan = y / x -> angle = arctan y / x "=" atan2(y, x)
    angle_rad = np.arctan2(a_y, a_x)
    angle_deg = np.rad2deg(angle_rad)

    atoms.rotate(-angle_deg, [0, 0, 1], rotate_cell=True)

    assert np.isclose(atoms.cell[0, 1], 0, atol=atol, rtol=0), atoms.cell
    assert np.isclose(atoms.cell[0, 2], 0, atol=atol, rtol=0), atoms.cell


def _align_b_onto_xy(atoms, atol):
    assert np.isclose(atoms.cell[0, 1], 0, atol=atol, rtol=0)  # make sure a is along x
    assert np.isclose(atoms.cell[0, 2], 0, atol=atol, rtol=0)  # make sure a is along x

    # rotate so that b is in xy plane
    # project b onto the yz-plane
    b = atoms.cell[1]
    b_y = b[1]
    b_z = b[2]
    angle_rad = np.arctan2(b_z, b_y)
    angle_deg = np.rad2deg(angle_rad)

    atoms.rotate(-angle_deg, [1, 0, 0], rotate_cell=True)

    assert np.isclose(atoms.cell[0, 1], 0, atol=atol, rtol=0)  # make sure a is in xy-plane
    assert np.isclose(atoms.cell[0, 2], 0, atol=atol, rtol=0)  # make sure a is in xy-plane
    assert np.isclose(atoms.cell[1, 2], 0, atol=atol, rtol=0), atoms.cell
