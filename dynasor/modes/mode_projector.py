from __future__ import annotations
import numpy as np
import warnings
import pickle
import itertools

from ase.calculators.singlepoint import SinglePointCalculator
from ase.units import fs
from ase import Atoms

from dynasor.units import radians_per_fs_to_THz
from .tools import get_dynamical_matrix, group_eigvals, symmetrize_eigenvectors, get_displacements
from .qpoint import QPoint
from .atoms import Prim, Supercell
from ..qpoints.tools import get_commensurate_lattice_points


class ModeProjector:
    """
    The :class:`ModeProjector` maps between real atomic displacements `u` and
    complex mode coordinates `Q`.

    Some special python methods are implemented. The `__str__` and `__repr__`
    provides useful info. :class:`QPoint` objects are representations of a
    single q-point and associated information and can be accessed either by
    call providing a reduced wavevector

    >>> mp((1/2, 0, 0))  # doctest: +SKIP

    or by index corresponding to reduced q-point accessible from
    :attr:`ModeProjector.q_reduced`

    >>> mp[2]  # doctest: +SKIP

    In addition to mode coordinates `Q` the class can also map the atomic
    velocities `v` to mode momenta `P` as well as atomic forces `f` to mode
    forces `F`. The class can also map back.  This mapping is done using
    getters and setters.  Internally only Q, P and F are stored

    >>> Q = mp.get_Q()  # doctest: +SKIP
    >>> mp.set_u(u)  # doctest: +SKIP

    In addition, the forces corresponding to the harmonic forces can be
    accessed by :meth:`ModeProjector.get_f_harmonic()` and
    :meth:`ModeProjector.get_F_harmonic()`. For ASE Atoms objects the
    displacments etc. can be updated and applied by

    >>> mp.update_from_atoms(atoms)  # doctest: +SKIP
    >>> atoms = mp.get_atoms(harmonic_forces=False)  # doctest: +SKIP

    The shapes for each property uses the follwoing varaibles

    * N: Number of atoms in supercell
    * Nu: unit cells (=N/Np)
    * Np: primitive basis atoms (=N/Nu)
    * Nb: bands (=Np*3)
    * Nq: q-points (=Nu)

    Please consult the documentation or the specific getters and setters to see
    the exact transformations used.

    Units
    -----
    The internal units in dynasor is Å, fs and eV. All frequencies are angular
    (a.k.a. "physicist's" convention with 2π included). These are the units
    dynasor will expect and return. In e.g. print functions conventional units
    like fs, Å, THz, Da, meV will often be used.

    Mass: The internal unit choice means that the mass unit is not Dalton
    (a.k.a. amu or u) but instead 0.009648533290731906 Da. We denote this with
    dmu (dynasor mass unit). In other words one Dalton is equal to
    103.64269572045423 dmu units.  This is typically not a problem since the
    masses provided via e.g.  ASE Atoms objects are converted internally.

    Waves: dynasor will always communicate and expect spatial (angular)
    frequencies in rad/Å and temporal (angular) frequencies in rad/fs. This
    follows the "physicist's" convention as the factor of 2π is included in the
    wave vectors. For instance the wavelength is given by λ=2π/q.

    Mode amplitudes: Usually the mode amplitudes is communicated in Å√dmu
    a.k.a. fs√eV.

    Velocities: For modes the momenta (Å√dmu/fs or just √eV) is used but for
    atoms the velocities (Å/fs) are used.

    Mode forces: The force is just defined as the rate of change of the
    canonical momenta so the unit would be Å√dmu/fs² (or √eV/fs)

    Internal arrays
    ---------------
    For the curious, the internal data arrays are

    * primitive, supercell, force_constants (input)
    * _q, q_minus (reduced q-points and which q-points are related by inversion)
    * _D, _w2, _W (dynamical matrices, frequencies (ev/Å²/Da), polarization vectors)
    * _X (eigenmodes which are mass weighted "polarization vectors" in the supercell)
    """
    def __init__(self, primitive, supercell, force_constants):
        """The mode projector is initialized by a primitive cell and a
        supercell as well as harmonic force constants.

        The force constants are assumed to be in units of eV/Å² as returned
        from phonopy. Be careful about the permutations when working with force
        constants and atoms object from different codes.

        Parameters
        ----------
        primitive: Atoms
            Primitive cell. Note that the masses are stored internally as
            Dalton in ASE but will be converted to dynasor's internal
            mass unit (dmu) internally in dynasor.
        supercell: Atoms
            Ideal supercell corresponding to the force constants.
        force_constants
            Force constants for the supercell in eV/Å² (N, N, 3, 3) wher N is
            `len(supercell)`
        """
        if len(primitive) == len(supercell):
            warnings.warn('Primitive and supercell have the same size')
        elif len(primitive) > len(supercell):
            raise ValueError('Primitive cell larger than supercell')
        elif not (len(supercell) / len(primitive)).is_integer():
            raise ValueError('supercell size is not multiple of primitive size')

        if len(supercell) != len(force_constants):
            raise ValueError('force constants shape is not compatible with supercell size')

        if force_constants.shape != (len(supercell), len(supercell), 3, 3):
            raise ValueError('force constants shape should be (N, N, 3, 3)')

        self.primitive = Prim(primitive)
        self.supercell = Supercell(supercell, primitive)
        self.force_constants = force_constants

        # Find q-points in reduced primitive cell coordinates
        q_integer = get_commensurate_lattice_points(self.supercell.P.T)
        q_reduced = np.dot(q_integer, self.supercell.P_inv.T)
        self._q = np.array(sorted(tuple(q) for q in q_reduced))

        # The equivalent q-point corresponding to -q
        self._q_minus = [[tuple(q) for q in self._q].index(tuple((-q) % 1)) for q in self._q]

        # Construct dynamical matrix and diagonalize at each q-point
        self._D, self._w2, self._W = [], [], []
        for qi, q in enumerate(self._q):

            D = get_dynamical_matrix(
                    self.force_constants, self.supercell.offsets, self.supercell.indices,
                    q.astype(np.float64))

            if qi == self.q_minus[qi]:
                assert np.allclose(D.imag, 0)
                D = D.real

            D = np.einsum('ijab,i,j->ijab',
                          D, self.primitive.masses**-0.5, self.primitive.masses**-0.5)
            D_matrix = D.transpose(0, 2, 1, 3).reshape(-1, self.primitive.n_atoms * 3)
            assert np.allclose(D_matrix, D_matrix.T.conj())
            w2, W = np.linalg.eigh(D_matrix)
            W = W.T.reshape(-1, self.primitive.n_atoms, 3)

            self._D.append(D)
            self._w2.append(w2)
            self._W.append(W)

        self._D = np.array(self._D)
        self._w2 = np.array(self._w2)
        self._W = np.array(self._W)

        # Post check basic symmetries, group eigenvalues and try to make degenerate modes nicer
        for q, q_minus in enumerate(self.q_minus):
            q_minus = self.q_minus[q]

            assert np.allclose(self._D[q], self._D[q_minus].conj())
            assert np.allclose(self._w2[q], self._w2[q_minus])

            # tolerances for grouping and sorting eigenvalues and eigenvectors
            round_decimals = 12
            tolerance = 10**(-round_decimals)

            for group in group_eigvals(self._w2[q], tolerance**0.5):
                W = symmetrize_eigenvectors(self._W[q, group])

                # Try to order them
                W_sort = W.copy().transpose(0, 2, 1).reshape(len(W), -1)
                # abs is because we want to consider the magnitude
                # - (minus) basically reverts the sort order to place largest first
                # T is just because how lexsort works, we want to consider each
                #     atom and direction as a key for the bands
                # -1 is because we want to make the x-direction of the first
                #     atom the mist significant key
                # At the end the first band should have the largest magnitude
                #     for the first atom in x
                argsort = np.lexsort(np.round(-np.abs(W_sort).T[::-1], round_decimals))
                self._W[q, group] = W[argsort]

            self._W[q_minus] = self._W[q].conj()

        # Construct supercell projection matrix
        # q_ks = X_ksna u_na
        self._X = np.zeros((len(self._q), self.primitive.n_atoms * 3, self.supercell.n_atoms, 3),
                           dtype=np.complex128)

        for index in range(self.supercell.n_atoms):
            i, N = self.supercell.indices[index], self.supercell.offsets[index]
            for q, s, a in itertools.product(
                    range(len(self._q)), range(self.primitive.n_atoms * 3), range(3)):
                phase = np.exp(-1j * 2*np.pi * self._q[q] @ N)
                self._X[q, s, index, a] = (
                        self.primitive.masses[i]**0.5 * phase * self._W[q, s, i, a].conj())
        self._X /= (self.supercell.n_atoms / self.primitive.n_atoms)**0.5

        # Init arrays to hold Q, P and F
        self._Q = np.zeros((len(self._q), self.primitive.n_atoms*3), dtype=np.complex128)
        self._P = np.zeros_like(self._Q)
        self._F = np.zeros_like(self._Q)

    # Dunder functions
    def __str__(self):
        strings = ['### ModeProjector ###']
        strings += [f'{self.supercell}']
        strings += [f'{self.primitive}']
        string = '\n'.join(strings)

        # ASCII DOS!
        width = 80
        height = 24
        dos = np.full((height, width), ' ')

        THz = self.omegas * radians_per_fs_to_THz

        hist, bins = np.histogram(THz.flat, bins=width)

        for i, h in enumerate(hist):
            dos[np.round(h * (height - 1) / hist.max()).astype(int), i] = '+'  # '·' or 'x'

        dos = dos[::-1]

        dos[-1, dos[-1] == ' '] = '-'

        dos = '\n'.join([''.join(d) for d in dos])

        string += f'\n{dos}'

        string += f'\n|{THz.min():<10.2f} THz' + ' '*(width - 26) + f'{THz.max():>10.2f}|'

        return string

    def __repr__(self):
        return str(self)

    def __getitem__(self, q) -> QPoint:
        """Returns the q-point object based on its index"""
        if q < 0 or q >= len(self._q):
            raise IndexError
        return QPoint(q, self)

    def __call__(self, qpoint) -> QPoint:
        """Tries to find a matching q-point based on reduced coordinate"""
        qpoint = np.array(qpoint).astype(np.float64) % 1
        for q, qpoint2 in enumerate(np.array(self._q).astype(np.float64)):
            if np.allclose(qpoint, qpoint2):
                return QPoint(q, self)
        raise ValueError('qpoint not compatible, check mp.q_reduced')

    # Getters ans setters for internal mutable arrays
    def get_Q(self):
        """The mode coordinate amplitudes in Å√dmu"""
        return self._Q.copy()

    def get_P(self):
        """The mode momentum amplitudes in √eV"""
        return self._P.copy()

    def get_F(self):
        """The mode force amplitudes in eV/Å√dmu"""
        return self._F.copy()

    def get_F_harmonic(self):
        """The harmonic mode forces, computed as -omega^2 * Q, in eV/Å√dmu"""
        return -self._w2 * self.get_Q()

    def set_Q(self, Q):
        """Sets the internal mode coordinates Q

        The functions ensures Q(-q)=Q*(q)
        """
        # This ensures that stuff like mp.set_Q(0) works while not updating the
        # array until the assert
        Q_new = self.get_Q()
        Q_new[:] = Q
        if not np.allclose(np.conjugate(Q_new), Q_new[self.q_minus]):
            raise ValueError('Supplied Q does not fulfill Q(-q) = Q(q)*')
        self._Q[:] = Q_new

    def set_P(self, P):
        """Sets the internal mode momenta P.

        The functions ensures :math:`P(-q)=P^*(q)`
        """
        P_new = self.get_P()
        P_new[:] = P
        if not np.allclose(np.conjugate(P_new), P_new[self.q_minus]):
            raise ValueError('Supplied P does not fulfill P(-q) = P(q)*')
        self._P[:] = P_new

    def set_F(self, F):
        """Sets the internal mode forces F.

        The functions ensures :math:`F(-q)=F^*(q)`
        """
        F_new = self.get_F()
        F_new[:] = F
        if not np.allclose(np.conjugate(F_new), F_new[self.q_minus]):
            raise ValueError('Supplied F does not fulfill F(-q) = F(q)*')
        self._F[:] = F_new

    def get_u(self):
        """The atomic displacements in Å"""
        u = np.einsum('ksna,ks,n->na', self._X.conj(), self._Q, 1 / self.supercell.masses)
        assert np.allclose(u.imag, 0)
        return u.real

    def get_v(self):
        """The atomic velocities in Å/fs"""
        v = np.einsum('ksna,ks,n->na', self._X, self._P, 1 / self.supercell.masses)
        assert np.allclose(v.imag, 0)
        return v.real

    def get_f(self):
        """The atomic forces in eV/Å"""
        f = np.einsum('ksna,ks->na', self._X, self._F)
        assert np.allclose(f.imag, 0)
        return f.real

    def get_f_harmonic(self):
        """The harmonic atomic forces for the current displacements"""
        F_harmonic = self.get_F_harmonic()
        f_harmonic = np.einsum('ksna,ks->na', self._X, F_harmonic)
        assert np.allclose(f_harmonic.imag, 0)
        return f_harmonic.real

    def set_u(self, u):
        """Sets the internal mode coordinates Q given the atomic displacements u

        `Q = X * u`

        Parameters
        ----------
        u
            The atomic displacements in Å
        """
        Q = np.einsum('ksna,na->ks', self._X, u)
        self.set_Q(Q)

    def set_v(self, v):
        """Sets the internal mode momenta P given the velocities v

        `P = X.conj() * v`

        Parameters
        ----------
        v
            The atomic velocities in fs/Å
        """
        P = np.einsum('ksna,na->ks', self._X.conj(), v)
        self.set_P(P)

    def set_f(self, f):
        """Sets the internal mode forces F given the forces f

        `F = X.conj() * f / m`

        Parameters
        ----------
        f
            The atomic forces in eV/Å
        """
        F = np.einsum('ksna,na,n->ks', self._X.conj(), f, 1 / self.supercell.masses)
        self.set_F(F)

    # Convenience functions to handle ASE Atoms objects
    def get_atoms(self, harmonic_forces=False) -> Atoms:
        """Returns ase Atoms object with displacement, velocities, forces and harmonic energies set

        Parameters
        ----------
        harmonic_forces
            Whether the forces should be taken from the internal `F` or via `-w2 Q`
        """
        atoms = self.supercell.to_ase()
        atoms.positions += self.get_u()
        atoms.set_velocities(self.get_v() / fs)
        E = self.potential_energies.sum()
        f = self.get_f_harmonic() if harmonic_forces else self.get_f()

        atoms.calc = SinglePointCalculator(
                energy=E, forces=f, stress=None, magmoms=None, atoms=atoms)

        return atoms

    def update_from_atoms(self, atoms):
        """Updates the ModeProjector with displacments, velocities and forces
        from an ASE Atoms object

        Checks for attached calculator in first place and next for forces array.

        If no data sets corresponding array to zeros.

        The masses and velocities are converted to dynasor units internally.
        """

        u = get_displacements(atoms, self.supercell)
        if np.max(np.abs(u)) > 2.0:
            warnings.warn('Displacements larger than 2Å. Is the atoms object permuted?')
        self.set_u(u)
        self.set_v(atoms.get_velocities() * fs)
        try:
            self.set_f(atoms.get_forces())
        except RuntimeError:
            if 'forces' in atoms.arrays:
                self.set_f(atoms.arrays['forces'])
            else:
                self.set_f(np.zeros_like(atoms.positions))

    # properties
    @property
    def q_minus(self):
        """The index of the corresponding counter-propagating mode (-q)"""
        return self._q_minus.copy()

    @property
    def q_reduced(self):
        """The q-points in reduced coordinates.

        A zone boundary mode would be e.g. (1/2, 0, 0)
        """
        return self._q.astype(float)

    @property
    def q_cartesian(self):
        """The q-points in cartesian coordinates with unit of rad/Å (2π included)"""
        return 2*np.pi * self.q_reduced @ self.primitive.inv_cell

    @property
    def omegas(self):
        """The frequencies of each mode in (rad/fs).

        Negative frequencies correspond to imaginary frequencies
        """
        return np.sign(self._w2) * np.sqrt(np.abs(self._w2))

    @property
    def polarizations(self):
        """The polarization vectors for each mode (Nq, Nb, Np, 3)"""
        return self._W

    @property
    def eigenmodes(self):
        """The eigenmodes in the supercell as (Nq, Nb, N, 3)-array

        The eigenmodes with masses included so that

        `Q = X u`

        where u is the supercell displacements
        """
        return self._X

    @property
    def potential_energies(self):
        """Potential energy per mode as (Nq, Nb)-array

        The potential energies are defined as :math:`1/2QQ^*` and should equal
        :math:`1/2k_BT` in equilibrium for a harmonic system.`
        """
        return 1/2 * np.abs(self._Q)**2 * self._w2

    @property
    def kinetic_energies(self):
        """Kinetic energy per mode as (Nq, Nb)-array

        The kinetic energies are defined as :math:`1/2PP^*`. Should equal
        :math:`1/2k_BT` in equilibrium.`
        """
        return 1/2 * np.abs(self._P)**2

    @property
    def virial_energies(self):
        """The virial energies per mode as (Nq, Nb)-array

        The virial energies are defined here as :math:`-1/2QF` which should have an
        expectation value of :math:`1/2k_BT` per mode in equilibrium. For a harmonic
        system this is simply equal to the potential energy. This means that
        that the virial energy can be used to monitor the anharmonicity or
        define a measure of the potential energy.
        """
        return -1/2 * self._Q * self._F

    def write(self, file_name: str):
        """Uses pickle to write mp to file"""
        with open(file_name, 'wb') as f:
            pickle.dump(self, f)

    @classmethod
    def read(cls, file_name: str) -> ModeProjector:
        """Return mp instance from pickle file that was saved using mp.write"""
        with open(file_name, 'rb') as f:
            mp = pickle.load(f)
        return mp
