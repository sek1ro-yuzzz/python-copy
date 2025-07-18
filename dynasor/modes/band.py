from __future__ import annotations
from typing import TYPE_CHECKING
from ..units import radians_per_fs_to_THz
from .tools import make_table
from .complex_coordinate import Q, P, F

if TYPE_CHECKING:
    from .qpoint import QPoint


class Band:
    """Represents properties of a single band belonging to a q-point

    To modify the coordinates `Q`, `P` and `F` use `band.Q` and see doc for
    :class:`dynasor.modes.complex_coordinate.ComplexCoordinate`
    """
    def __init__(self, q_index, band_index, mp):
        self._q = q_index
        self._s = band_index
        self._mp = mp

    def __repr__(self):
        return str(self)

    def __str__(self):
        string = f"""### Band ###
Band index:     {self.index}
q-point:        [{self.qpoint.index}] {self.qpoint.q_reduced}
Frequency:      {self.omega * radians_per_fs_to_THz:.2f} THz
Polarization:
"""
        string += make_table(self.polarization)
        string += f"""
              {'Q':<20}{'P':<20}{'F':<20}
Value:        {self.Q.complex:<20.2f}{self.P.complex:<20.2f}{self.F.complex:<20.2f}
Amplitude:    {self.Q.amplitude:<20.2f}{self.P.amplitude:<20.2f}{self.F.amplitude:<20.2f}
Phase:        {self.Q.phase:<20.2f}{self.P.phase:<20.2f}{self.F.phase:<20.2f}
Energy:       {self.potential_energy:<20.2f}{self.kinetic_energy:<20.2f}{self.virial_energy:<20.2f}
"""
        return string

    @property
    def index(self):
        """The index of the band at the q-point"""
        return self._s

    @property
    def qpoint(self) -> 'QPoint':
        """The q-point to which the band belongs"""
        return self._mp[self._q]

    @property
    def omega(self):
        """Slice, see doc for ModeProjector"""
        return self.qpoint.omegas[self.index]

    @property
    def polarization(self):
        """Slice, see doc for ModeProjector"""
        return self.qpoint.polarizations[self.index]

    @property
    def eigenmode(self):
        """Slice, see doc for ModeProjector"""
        return self.qpoint.eigenmodes[self.index]

    @property
    def potential_energy(self):
        """Slice, see doc for ModeProjector"""
        return self.qpoint.potential_energies[self.index]

    @property
    def kinetic_energy(self):
        """Slice, see doc for ModeProjector"""
        return self.qpoint.kinetic_energies[self.index]

    @property
    def virial_energy(self):
        """Slice, see doc for ModeProjector"""
        return self.qpoint.virial_energies[self.index]

    @property
    def Q(self) -> Q:
        """The mode coordinate"""
        return Q(self.qpoint.index, self.index, self._mp)

    @property
    def P(self) -> P:
        """The mode momentum"""
        return P(self.qpoint.index, self.index, self._mp)

    @property
    def F(self) -> F:
        """The mode force"""
        return F(self.qpoint.index, self.index, self._mp)
