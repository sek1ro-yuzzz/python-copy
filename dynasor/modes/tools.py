import numbers
import numba
from ase.geometry import get_distances
import numpy as np
from fractions import Fraction
from ase.geometry import find_mic


# Geometry
def get_displacements(supercell, ideal, check_mic=True):
    """Gets displacement from a supercell that might have wrapped positions"""
    u = supercell.positions - ideal.positions
    return get_wrapped_displacements(u, ideal.cell, check_mic=True)


def get_wrapped_displacements(u, cell, check_mic=True):
    """wraps displacements using mic"""
    if check_mic:
        u, _ = find_mic(u, cell)
    return u


def find_permutation(atoms, atoms_ref):
    """ Returns the best permutation of atoms for mapping one
    configuration onto another.

    Parameters
    ----------
    atoms
        configuration to be permuted
    atoms_ref
        configuration onto which to map

    Examples
    --------
    After obtaining the permutation via ``p = find_permutation(atoms1, atoms2)``
    the reordered structure ``atoms1[p]`` will give the closest match
    to ``atoms2``.
    """
    assert np.linalg.norm(atoms.cell - atoms_ref.cell) < 1e-6
    permutation = []
    for i in range(len(atoms_ref)):
        dist_row = get_distances(
            atoms.positions, atoms_ref.positions[i], cell=atoms_ref.cell, pbc=True)[1][:, 0]
        permutation.append(np.argmin(dist_row))

    if len(set(permutation)) != len(permutation):
        raise Exception('Duplicates in permutation')
    for i, p in enumerate(permutation):
        if atoms[p].symbol != atoms_ref[i].symbol:
            raise Exception('Matching lattice sites have different occupation')
    return permutation


# Linalg
def trace(A):
    return A[0, 0] + A[1, 1] + A[2, 2]


def det(A):
    if len(A) == 2:
        return A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]
    d = 0
    for i, B in enumerate(A[0]):
        minor = np.hstack([A[1:, :i], A[1:, i+1:]])
        d += (-1)**i * B * det(minor)
    return d


def inv(A, as_fraction=True):
    """
    Inverts A matrix

    Parameters
    ----------
    A
        array to be inverted
    as_fraction
        boolean which determines if inverted matrix is returned as Fraction matrix or real matrix
    """

    detx2 = det(A) * 2

    rest = (trace(A)**2 - trace(A @ A)) * np.diag([1, 1, 1]) - 2 * A * trace(A) + 2 * A @ A

    if detx2 < 0:
        detx2 = -detx2
        rest = -rest

    A_inv = np.array([Fraction(n, detx2) for n in rest.flat]).reshape(3, 3)

    assert np.all(A @ A_inv == np.eye(len(A)))

    # return either real matrix or as fraction matrix
    A_inv_real = np.reshape([float(x) for x in A_inv.flat], A_inv.shape)
    assert np.allclose(A_inv_real, np.linalg.inv(A))

    if as_fraction:
        return A_inv
    else:
        return A_inv_real


def symmetrize_eigenvectors(eigenvectors, cell=None, max_iter=1000, method='varimax'):
    """Takes a set of vectors and tries to make them nice

    Parameters
    ----------
    eigenvectors
        If there are n degenerate eigenvectors and m atoms in the basis the
        array should be (n,m,3)
    cell
        If default None nothing is done to the cartesian directions but a cell
        can be provided so the directions are in scaled coordinates instead.
    max_iter
        The number of iterations in the symmetrization procedure
    method
        Can be 'varimax' or 'quartimax' or a parameter between 0: qartimax and
        1: varimax. Depending on the choice one gets e.g. Equamax, Parsimax, etc.

    """

    if cell is None:
        cell = np.eye(3)

    # s = band, i = basis, a = axis
    eigenvectors = np.einsum('sia,ab->sib', eigenvectors, np.linalg.inv(cell))

    components = eigenvectors.reshape(len(eigenvectors), -1).T

    rotation_matrix = factor_analysis(components, iterations=max_iter, method=method)

    new_eigenvectors = np.dot(components, rotation_matrix).T

    new_eigenvectors = new_eigenvectors.reshape(len(new_eigenvectors), -1, 3)

    new_eigenvectors = np.einsum('sia,ab->sib', new_eigenvectors, cell)

    return new_eigenvectors


def factor_analysis(L, iterations=1000, method='varimax'):
    """Performs factor analysis on L finding rotation matrix R such that L @ R = L' is simple

    In the future consider using the scikit learn methods directly but beware
    the changes need to accomodate complex numbers

    References:
        * Sparse Modeling of Landmark and Texture Variability using the Orthomax Criterion
        Mikkel B. Stegmann, Karl Sj¨ostrand, Rasmus Larsen
        http://www2.imm.dtu.dk/pubdb/edoc/imm4041.pdf

        * http://www.cs.ucl.ac.uk/staff/d.barber/brml,

        * https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.FactorAnalysis.html   # noqa

        * https://stats.stackexchange.com/questions/185216/factor-rotation-methods-varimax-quartimax-oblimin-etc-what-do-the-names  # noqa

    """

    nrow, ncol = L.shape
    R = np.eye(ncol)

    if method == 'varimax':
        gamma = 1
    elif method == 'quartimax':
        gamma = 0
    else:
        gamma = method

    for _ in range(iterations):
        LR = np.dot(L, R)
        grad = LR * (np.abs(LR)**2 - gamma * np.mean(np.abs(LR)**2, axis=0))
        G = L.T.conj() @ grad
        u, s, vh = np.linalg.svd(G)
        R = u @ vh

    return R


def group_eigvals(vals, tol=1e-6):
    assert sorted(vals) == list(vals), vals

    groups = [[0]]
    for i in range(1, len(vals)):
        if np.abs(vals[i] - vals[i-1]) < tol:
            groups[-1].append(i)
        else:
            groups.append([i])
    return groups


# misc
def as_fraction(not_fraction):

    if isinstance(not_fraction, numbers.Number):
        return Fraction(not_fraction)

    if isinstance(not_fraction, np.ndarray):
        arr = np.array([Fraction(n) for n in not_fraction.flat])
        return arr.reshape(not_fraction.shape)

    if isinstance(not_fraction, tuple):
        arr = (Fraction(n) for n in not_fraction)
        return arr

    if isinstance(not_fraction, list):
        arr = [Fraction(n) for n in not_fraction]
        return arr


@numba.njit
def get_dynamical_matrix(fc, offsets, indices, q):

    n = indices.max() + 1
    N = len(fc)
    D = np.zeros(shape=(n, n, 3, 3), dtype=np.complex128)
    for ia in range(n):
        for a in range(N):
            if ia != indices[a]:
                continue
            na = offsets[a]
            for b in range(N):
                ib = indices[b]
                nb = offsets[b]
                dn = (nb - na).astype(np.float64)
                D[ia, ib] += fc[a, b] * np.exp(2j*np.pi * np.dot(dn, q))
            break
    return D


# For debug, this is a slower but perhaps more accurate variant, if the fc
# obeys translational invariance this should give the same result
# @numba.njit
# def get_dynamical_matrix_full(fc, offsets, indices, q):
#
#     n = indices.max() + 1
#     N = len(fc)
#     D = np.zeros(shape=(n, n, 3, 3), dtype=np.complex128)
#     for I in range(N):
#         i = indices[I]
#         m = offsets[I]
#         for J in range(N):
#             j = indices[J]
#             n = offsets[J]
#
#             off = (n - m).astype(np.float64)
#
#             phase = np.exp(1j * 2*np.pi * np.dot(off, q))
#
#             D[i, j] += fc[I, J] * phase
#
#     D /= (N / D.shape[0])
#
#     return D


def make_table(M):
    rows = []
    for r in M:
        rows.append(''.join(f'{e:<20.2f}' for e in r))
    return '\n'.join(rows)
