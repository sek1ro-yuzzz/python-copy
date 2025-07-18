"""
The mode projection functionality in dynasor is mainly handled by two objects:
the :class:`dynasor.ModeProjector` class and :class:`dynasor.project_modes`
function.  From the :class:`dynasor.ModeProjector` there is access to
data-objects representing a q-point :class:`dynasor.modes.qpoint.QPoint` and
from the q-point there is access to an object representing a particular band
:class:`dynasor.modes.band.Band` at that q-point.  In addition, simple wrappers
around the coordinates `Q`, `P` and `F` exists via
:class:`dynasor.modes.complex_coordinate.ComplexCoordinate` to easily set the
amplitude and phase of a mode while preserving the :math:`Q(-q)=Q^*(q)`
symmetries.  Internally dynasor wraps the primitive cell
:class:`dynasor.modes.atoms.Prim` and supercell
:class:`dynasor.modes.atoms.Supercell`.  As a user only the
:class:`dynasor.ModeProjector` and :class:`dynasor.project_modes` should need
to be imported.
"""
