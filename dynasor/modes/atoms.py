import numpy as np

from ase.units import fs
from ase import Atoms

from .tools import inv
from ..units import Dalton_to_dmu


class DynasorAtoms:
    """Dynasor's representation of a structure"""
    def __init__(self, atoms: Atoms):
        """Initialized using an ASE Atoms object"""
        self._atoms = atoms

    @property
    def pos(self):
        """Cartesian positions"""
        return self._atoms.positions.copy()

    @property
    def positions(self):
        """Cartesian positions"""
        return self.pos

    @property
    def spos(self):
        """Reduced (or scaled) positions of atoms"""
        return self._atoms.get_scaled_positions()

    @property
    def scaled_positions(self):
        """Reduced (or scaled) positions of atoms"""
        return self.spos

    @property
    def cell(self):
        """Cell of atoms with cell vectors as rows"""
        return self._atoms.cell.array.copy()

    @property
    def inv_cell(self):
        """The inverse cell transpose so the inverse cell vectors are rows, no 2pi"""
        return np.linalg.inv(self._atoms.cell.array).T

    @property
    def numbers(self):
        """Chemical number for each atom, e.g. 1 for H, 2 for He etc."""
        return self._atoms.numbers.copy()

    @property
    def masses(self):
        """Masses of atoms in dmu"""
        return self._atoms.get_masses() / fs**2  # In eVfs²/Å²

    @property
    def volume(self):
        """Volume of cell"""
        return self._atoms.get_volume()

    @property
    def n_atoms(self):
        """Number of atoms"""
        return len(self._atoms)

    @property
    def symbols(self):
        """List of chemical symbol for each element"""
        return list(self._atoms.symbols)

    def to_ase(self):
        """Converts the internal Atoms to ASE Atoms"""
        return Atoms(cell=self.cell, numbers=self.numbers, positions=self.positions, pbc=True)

    def __repr__(self):
        return str(self)


class Prim(DynasorAtoms):
    def __str__(self):
        strings = [f"""Primitive cell:
Number of atoms:        {self.n_atoms}
Volume:                 {self.volume:.3f}
Atomic species present: {set(self.symbols)}
Atomic numbers present: {set(self.numbers)}
Cell:
[[{self.cell[0, 0]:<20}, {self.cell[0, 1]:<20}, {self.cell[0, 2]:<20}],
 [{self.cell[1, 0]:<20}, {self.cell[1, 1]:<20}, {self.cell[1, 2]:<20}],
 [{self.cell[2, 0]:<20}, {self.cell[2, 1]:<20}, {self.cell[2, 2]:<20}]]
"""]
        strings.append(f"{'Ind':<5}{'Sym':<5}{'Num':<5}{'Mass (Da)':<10}{'x':<10}{'y':<10}{'z':<10}"
                       f"{'a':<10}{'b':<10}{'c':<10}")
        atom_s = []
        for i, p, sp, m, n, s in zip(
                range(self.n_atoms), self.positions, self.spos, self.masses / Dalton_to_dmu,
                self.numbers, [a.symbol for a in self.to_ase()]):
            atom_s.append(f'{i:<5}{s:<5}{n:<5}{m:<10.2f}{p[0]:<10.3f}{p[1]:<10.3f}{p[2]:<10.3f}'
                          f'{sp[0]:<10.3f}{sp[1]:<10.3f}{sp[2]:<10.3f}')

        strings = strings + atom_s

        string = '\n'.join(strings)

        return string


class Supercell(DynasorAtoms):
    """The supercell takes care of some mappings between the primitive and repeated structure

    In particular the P-matrix connecting the cells as well as the offset-index of each atom is
    calculated.

    Note that the positions can not be revovered as offset x cell + basis since the atoms gets
    wrapped.

    Parameters
    ----------
    supercell
        some ideal repetition of the prim cell and possible wrapping  as either ASEAtoms
    prim
        primitive cell as either ASEAtoms

    """

    def __init__(self, supercell, prim):
        self.prim = Prim(prim.copy())
        super().__init__(supercell)

        # determine P-matrix relating supercell to primitive cell
        from dynasor.tools.structures import get_P_matrix
        self._P = get_P_matrix(self.prim.cell, self.cell)  # P C = S
        self._P_inv = inv(self.P)

        # find the index and offsets for supercell using primitive as base unit
        from dynasor.tools.structures import get_offset_index
        self._offsets, self._indices = get_offset_index(prim, supercell, wrap=True)

    @property
    def P(self):
        """P-matrix is defined as dot(P, prim.cell) = supercell.cell"""
        return self._P.copy()

    @property
    def P_inv(self):
        """Inverse of P"""
        return self._P_inv.copy()

    @property
    def offsets(self):
        """The offset of each atom"""
        return self._offsets.copy()

    @property
    def indices(self):
        """The basis index of each atom"""
        return self._indices.copy()

    @property
    def n_cells(self):
        """Number of unit cells"""
        return self.n_atoms // self.prim.n_atoms

    def __str__(self):

        string = f"""Supercell:
Number of atoms:      {self.n_atoms}
Volume:               {self.volume:.3f}
Number of unit cells: {self.n_cells}
Cell:
[[{self.cell[0, 0]:<20}, {self.cell[0, 1]:<20}, {self.cell[0, 2]:<20}],
 [{self.cell[1, 0]:<20}, {self.cell[1, 1]:<20}, {self.cell[1, 2]:<20}],
 [{self.cell[2, 0]:<20}, {self.cell[2, 1]:<20}, {self.cell[2, 2]:<20}]]
P-matrix:
[[{self.P[0, 0]:<20}, {self.P[0, 1]:<20}, {self.P[0, 2]:<20}],
 [{self.P[1, 0]:<20}, {self.P[1, 1]:<20}, {self.P[1, 2]:<20}],
 [{self.P[2, 0]:<20}, {self.P[2, 1]:<20}, {self.P[2, 2]:<20}]]
{self.prim}
"""
        return string
