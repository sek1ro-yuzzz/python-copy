from __future__ import annotations
import numpy as np

from ..units import radians_per_fs_to_THz
from .band import Band


class QPoint:
    """Representation of a single q-point and properties

    The bands can be accessed by e.g. `qp[2]` to get band index 2 in the form
    of a :class:`dynasor.modes.band.Band` object
    """
    def __init__(self, q, mp):
        self._mp = mp
        self._q = q

    # Dunders

    def __str__(self):
        string = f"""### q-point ###
Index:       {self.index}
Reduced:     {self.q_reduced}
Cartesian:   [{self.q_cartesian[0]:.2f}, {self.q_cartesian[1]:.2f}, {self.q_cartesian[2]:.2f}] rad/Å
Wavenumber:  {self.wavenumber:.2f} rad/Å
Wavelength:  {self.wavelength:.2f} Å
q-minus:     {self.q_minus.index} {self.q_minus.q_reduced}
Real:        {self.is_real}"""
        string += ("""
Omegas:      [""" + ', '.join(f'{t * radians_per_fs_to_THz:.2f}' for t in self.omegas) + '] THz')
        return string

    def __repr__(self):
        return str(self)

    def __getitem__(self, s):
        if s >= self._mp.primitive.n_atoms * 3:
            raise IndexError
        return Band(self.index, s, self._mp)

    @property
    def q_minus(self) -> QPoint:
        """The corresponding counter-propagating mode"""
        return self._mp[self._mp.q_minus[self.index]]

    @property
    def polarizations(self):
        """Slice, see :class:`dynasor.ModeProjector`"""
        return self._mp.polarizations[self.index]

    @property
    def omegas(self):
        """Slice, see :class:`dynasor.ModeProjector`"""
        return self._mp.omegas[self.index]

    @property
    def is_real(self):
        """If the q-point has purely real mode coordinates, 'q=-q'"""
        return self.index == self.q_minus.index

    @property
    def index(self):
        """q-point index corresponding to :class:`dynasor.ModeProjector.q_reduced`"""
        return self._q

    @property
    def wavenumber(self):
        """Wavenumber of mode in rad/Å"""
        return np.linalg.norm(self._mp.q_cartesian[self.index])

    @property
    def wavelength(self):
        """Wavelength of mode in Å"""
        return np.inf if self.wavenumber == 0 else 2*np.pi / self.wavenumber

    @property
    def q_reduced(self):
        """Slice, see :class:`dynasor.ModeProjector`"""
        return self._mp.q_reduced[self.index]

    @property
    def q_cartesian(self):
        """Slice, see :class:`dynasor.ModeProjector`"""
        return self._mp.q_cartesian[self.index]

    @property
    def eigenmodes(self):
        """Slice, see :class:`dynasor.ModeProjector`"""
        return self._mp.eigenmodes[self.index]

    @property
    def potential_energies(self):
        """Slice, see :class:`dynasor.ModeProjector`"""
        return self._mp.potential_energies[self.index]

    @property
    def kinetic_energies(self):
        """Slice, see :class:`dynasor.ModeProjector`"""
        return self._mp.kinetic_energies[self.index]

    @property
    def virial_energies(self):
        """Slice, see :class:`dynasor.ModeProjector`"""
        return self._mp.virial_energies[self.index]
