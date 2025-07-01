"""This module contains all modifiers available in OVITO. See the introduction to learn more
about modifiers and their role in the data pipeline system. The following table lists the Python names of all modifier types that can be instantiated.
Please consult the OVITO user manual for a more in-depth description of what
each of these modifiers does.

============================================== =========================================
Python class name                              User interface name
============================================== =========================================
:py:class:`AcklandJonesModifier`               :guilabel:`Ackland-Jones analysis`
:py:class:`AffineTransformationModifier`       :guilabel:`Affine transformation`
:py:class:`AmbientOcclusionModifier`           :guilabel:`Ambient occlusion`
:py:class:`AssignColorModifier`                :guilabel:`Assign color`
:py:class:`AtomicStrainModifier`               :guilabel:`Atomic strain`
:py:class:`BondAnalysisModifier`               :guilabel:`Bond analysis`
:py:class:`CalculateDisplacementsModifier`     :guilabel:`Displacement vectors`
:py:class:`CentroSymmetryModifier`             :guilabel:`Centrosymmetry parameter`
:py:class:`ChillPlusModifier`                  :guilabel:`Chill+`
:py:class:`ClearSelectionModifier`             :guilabel:`Clear selection`
:py:class:`ClusterAnalysisModifier`            :guilabel:`Cluster analysis`
:py:class:`ColorByTypeModifier`                :guilabel:`Color by type`
:py:class:`ColorCodingModifier`                :guilabel:`Color coding`
:py:class:`CombineDatasetsModifier`            :guilabel:`Combine datasets`
:py:class:`CommonNeighborAnalysisModifier`     :guilabel:`Common neighbor analysis`
:py:class:`ComputePropertyModifier`            :guilabel:`Compute property`
:py:class:`ConstructSurfaceModifier`           :guilabel:`Construct surface mesh`
:py:class:`CoordinationAnalysisModifier`       :guilabel:`Coordination analysis`
:py:class:`CoordinationPolyhedraModifier`      :guilabel:`Coordination polyhedra`
:py:class:`CreateBondsModifier`                :guilabel:`Create bonds`
:py:class:`CreateIsosurfaceModifier`           :guilabel:`Create isosurface`
:py:class:`DeleteSelectedModifier`             :guilabel:`Delete selected`
:py:class:`DislocationAnalysisModifier`        :guilabel:`Dislocation analysis (DXA)`
:py:class:`ElasticStrainModifier`              :guilabel:`Elastic strain calculation`
:py:class:`ExpandSelectionModifier`            :guilabel:`Expand selection`
:py:class:`ExpressionSelectionModifier`        :guilabel:`Expression selection`
:py:class:`FreezePropertyModifier`             :guilabel:`Freeze property`
:py:class:`GenerateTrajectoryLinesModifier`    :guilabel:`Generate trajectory lines`
:py:class:`GrainSegmentationModifier`          :guilabel:`Grain segmentation`
:py:class:`HistogramModifier`                  :guilabel:`Histogram`
:py:class:`IdentifyFCCPlanarFaultsModifier`    :guilabel:`Identify fcc planar faults`
:py:class:`IdentifyDiamondModifier`            :guilabel:`Identify diamond structure`
:py:class:`InvertSelectionModifier`            :guilabel:`Invert selection`
:py:class:`LoadTrajectoryModifier`             :guilabel:`Load trajectory`
:py:class:`PolyhedralTemplateMatchingModifier` :guilabel:`Polyhedral template matching`
:py:class:`PythonModifier`                     :guilabel:`Python script`
:py:class:`RenderLAMMPSRegionsModifier`        :guilabel:`Render LAMMPS regions`
:py:class:`ReplicateModifier`                  :guilabel:`Replicate`
:py:class:`SelectTypeModifier`                 :guilabel:`Select type`
:py:class:`SliceModifier`                      :guilabel:`Slice`
:py:class:`SmoothTrajectoryModifier`           :guilabel:`Smooth trajectory`
:py:class:`SpatialBinningModifier`             :guilabel:`Spatial binning`
:py:class:`SpatialCorrelationFunctionModifier` :guilabel:`Spatial correlation function`
:py:class:`TimeAveragingModifier`              :guilabel:`Time averaging`
:py:class:`TimeSeriesModifier`                 :guilabel:`Time series`
:py:class:`UnwrapTrajectoriesModifier`         :guilabel:`Unwrap trajectories`
:py:class:`VoronoiAnalysisModifier`            :guilabel:`Voronoi analysis`
:py:class:`VoroTopModifier`                    :guilabel:`VoroTop analysis`
:py:class:`WignerSeitzAnalysisModifier`        :guilabel:`Wigner-Seitz defect analysis`
:py:class:`WrapPeriodicImagesModifier`         :guilabel:`Wrap at periodic boundaries`
============================================== =========================================

*Note that some modifiers of the graphical version of OVITO are missing from this list and are not accessible from Python scripts.
That is because they perform simple operations that can be accomplished equally well or even easier using other means in Python.*"""
__all__ = ['PythonModifier', 'PythonScriptModifier', 'TimeAveragingModifier', 'TimeSeriesModifier', 'SliceModifier', 'AffineTransformationModifier', 'ClearSelectionModifier', 'InvertSelectionModifier', 'ColorCodingModifier', 'AssignColorModifier', 'DeleteSelectedModifier', 'SelectTypeModifier', 'HistogramModifier', 'ScatterPlotModifier', 'ReplicateModifier', 'ExpressionSelectionModifier', 'FreezePropertyModifier', 'ManualSelectionModifier', 'ComputePropertyModifier', 'CombineDatasetsModifier', 'ColorByTypeModifier', 'CreateIsosurfaceModifier', 'AmbientOcclusionModifier', 'WrapPeriodicImagesModifier', 'ExpandSelectionModifier', 'StructureIdentificationModifier', 'CommonNeighborAnalysisModifier', 'AcklandJonesModifier', 'CreateBondsModifier', 'CentroSymmetryModifier', 'ClusterAnalysisModifier', 'CoordinationAnalysisModifier', 'CalculateDisplacementsModifier', 'AtomicStrainModifier', 'WignerSeitzAnalysisModifier', 'VoronoiAnalysisModifier', 'IdentifyDiamondModifier', 'LoadTrajectoryModifier', 'PolyhedralTemplateMatchingModifier', 'CoordinationPolyhedraModifier', 'SmoothTrajectoryModifier', 'GenerateTrajectoryLinesModifier', 'UnwrapTrajectoriesModifier', 'ChillPlusModifier', 'ConstructSurfaceModifier', 'CoordinationNumberModifier', 'InterpolateTrajectoryModifier', 'SpatialBinningModifier', 'BondAnalysisModifier', 'SpatialCorrelationFunctionModifier', 'DislocationAnalysisModifier', 'ElasticStrainModifier', 'GrainSegmentationModifier', 'IdentifyFCCPlanarFaultsModifier', 'RenderLAMMPSRegionsModifier', 'CalculateLocalEntropyFunction', 'ShrinkWrapSimulationBoxFunction', 'VoroTopModifier']
from __future__ import annotations
from typing import Tuple, Optional, Any, Union, MutableSet, Sequence, Callable, Generator, Mapping
import ovito.vis
import ovito.pipeline
import ovito.data
from ovito import ArrayLike
from numpy.typing import NDArray
import enum
import math
from dataclasses import dataclass

@dataclass(kw_only=True)
class StructureIdentificationModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

Abstract base class for all modifiers in OVITO that analyze the local neighborhood of particles to 
identify structural motives or crystalline structures. It contains parameter fields that are
common to all these modifiers:

  * :py:class:`AcklandJonesModifier`
  * :py:class:`ChillPlusModifier`
  * :py:class:`CommonNeighborAnalysisModifier`
  * :py:class:`IdentifyDiamondModifier`
  * :py:class:`PolyhedralTemplateMatchingModifier`
  * :py:class:`VoroTopModifier`"""
    color_by_type: bool = True
    'color_by_type() -> bool\n\nControls the coloring of particles by the modifier to indicate their identified structure type. See the documentation of the :py:attr:`structures` list on how to change the colors representing the different structure types recognized by the modifier. \n\nDefault: ``True``'
    only_selected: bool = False
    'only_selected() -> bool\n\nSet this to ``True`` to perform the analysis on selected particles only. Particles that are *not* selected will be treated as if they did not exist and are assigned to the "OTHER" structure category. Use a :py:class:`SelectTypeModifier` in your pipeline, for example, to restrict the structure identification to a sub-lattice formed by one species of particles in a multi-component system. \n\nDefault: ``False``'

    @property
    def structures(self) -> Sequence[ovito.data.ElementType]:
        """A list of :py:class:`ElementType` instances managed by this modifier, one for each recognized structural type. You can modify the type objects in this list to adjust the coloring and to turn the identification of certain structural types on or off. The ordering of the :py:class:`ParticleType` objects in this list corresponds to the numeric type IDs defined by the concrete structure identification modifiers. In the following code snippets, the :py:class:`CommonNeighborAnalysisModifier` serves as an example. 

Your can change the color of a type by setting its :py:attr:`color` property to a new RGB value:: 

   modifier = CommonNeighborAnalysisModifier()
   modifier.structures[CommonNeighborAnalysisModifier.Type.FCC].color = (0.2, 1.0, 0.8)
   modifier.structures[CommonNeighborAnalysisModifier.Type.HCP].color = (0.0, 0.4, 1.0)


To turn the identification of a particular structure type on or off, you set its :py:attr:`enabled` property:: 

   modifier.structures[CommonNeighborAnalysisModifier.Type.BCC].enabled = False
   modifier.structures[CommonNeighborAnalysisModifier.Type.ICO].enabled = True




."""
        ...

@dataclass(kw_only=True)
class AcklandJonesModifier(StructureIdentificationModifier):
    """Base: :py:class:`ovito.modifiers.StructureIdentificationModifier`

This modifier analyzes the local neighborhood of each particle to identify simple crystalline structures.
The structure identification is performed using the bond-angle classification method proposed by Ackland and Jones.
See the corresponding user manual page
for more information on this modifier.

Note that this class inherits several important parameter fields from its :py:class:`StructureIdentificationModifier`
base class.

Modifier inputs:

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Particle properties
      -
    * - ``Position``
      - The input coordinates of the particles.
    * - ``Selection``
      - The selection state of the input particles. Only needed if :py:attr:`~StructureIdentificationModifier.only_selected` is ``True``.

Modifier outputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Structure Type``
      - The structure type computed by the algorithm for each particle, encoded as an integer value:

        ============= =========================================================
        Numeric id    Python constant
        ============= =========================================================
        0             ``AcklandJonesModifier.Type.OTHER``
        1             ``AcklandJonesModifier.Type.FCC``
        2             ``AcklandJonesModifier.Type.HCP``
        3             ``AcklandJonesModifier.Type.BCC``
        4             ``AcklandJonesModifier.Type.ICO``
        ============= =========================================================
    * - ``Color``
      - Particle coloring to indicate the identified structure type for each particle; only if :py:attr:`~StructureIdentificationModifier.color_by_type` is ``True``.
        See the :py:attr:`~StructureIdentificationModifier.structures` array on how to customize the colors.

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Global attributes
      -
    * - ``AcklandJones.counts.OTHER``
      - Number of particles not matching any of the recognized structure types.
    * - ``AcklandJones.counts.FCC``
      - Number of particles identified as face-centered cubic.
    * - ``AcklandJones.counts.HCP``
      - Number of particles identified as hexagonal close packed.
    * - ``AcklandJones.counts.BCC``
      - Number of particles identified as body-centered cubic.
    * - ``AcklandJones.counts.ICO``
      - Number of particles identified as icosahedral.

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Data tables
      -
    * - ``structures``
      - A bar chart with the particle counts for each structure type identified by the modifier.
        You can retrieve this :py:class:`DataTable` from the `DataCollection.tables` dictionary."""

    class Type(enum.IntEnum):
        """"""
        OTHER = enum.auto()
        FCC = enum.auto()
        HCP = enum.auto()
        BCC = enum.auto()
        ICO = enum.auto()

@dataclass(kw_only=True)
class AffineTransformationModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

This modifier applies an affine transformation to data elements in order to move, rotate, shear or scale them.
See also the corresponding user manual page for more information.

Inputs:

The modifier can operate on any combination of the following data elements:

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Data element specifier
      - Description
    * - ``particles``
      - Transforms the ``Position`` particle property.
    * - ``vector_properties``
      - Transforms all vectorial properties having 3 components and an attached :py:class:`VectorVis` element, e.g. the particle properties ``Velocity``, ``Force``, and ``Displacement``.
    * - ``cell``
      - Transforms the :py:class:`SimulationCell`.
    * - ``lines``
      - Transforms the vertices of all :py:class:`Lines` objects.
    * - ``surfaces``
      - Transforms the vertices of all :py:class:`SurfaceMesh` and :py:class:`TriangleMesh` objects.
    * - ``dislocations``
      - Transforms all dislocation lines in a :py:class:`DislocationNetwork`.
    * - ``vectors``
      - Transforms the base points and directions of all :py:class:`Vectors` objects.
    * - ``voxels``
      - Transforms the spatial :py:attr:`domain` of all :py:class:`VoxelGrid` objects.

The modifier will act on all of these by default. You can restrict the transform to a subset of data objects by setting the :py:attr:`operate_on` field.

Examples:

The following code applies a simple shear transformation to all particle coordinates. The shear transformation is specified
as a 3x3 matrix with ones on the matrix diagonal and an off-diagonal element that is non-zero.
No translation is applied to the particles. Thus all elements in the fourth column of the extended 3x4 affine transformation
matrix are set to zero:

```python
  xy_shear = 0.05
  mod = AffineTransformationModifier(
          operate_on = {'particles'}, # Transform particles but not the box.
          transformation = [[1, xy_shear, 0, 0],
                            [0,        1, 0, 0],
                            [0,        0, 1, 0]])
```

Note that the modifier itself only supports static transformations, which
remain constant over the entire trajectory. However, is possible to employ the :py:class:`AffineTransformationModifier`
in a user-defined modifier function, which calculates the transformation matrix dynamically
at each animation frame:

```python
  import numpy as np
  
  def rotate(frame, data):
      theta = np.deg2rad(frame * 5.0)  # time-dependent angle of rotation
      tm = [[np.cos(theta), -np.sin(theta), 0, 0],
            [np.sin(theta),  np.cos(theta), 0, 0],
            [            0,              0, 1, 0]]
      # Execute AffineTransformationModifier as a sub-operation:
      data.apply(AffineTransformationModifier(transformation = tm))
  
  pipeline.modifiers.append(rotate)
```"""
    only_selected: bool = False
    'only_selected() -> bool\n\nControls whether the modifier should transform only the subset of currently selected elements (e.g. selected particles). For data object types that do not support selections, e.g. the simulation cell, this option has no effect. \n\nDefault: ``False``'
    operate_on: MutableSet[str] = {'particles', 'vector_properties', 'cell', 'surfaces', 'dislocations', 'lines', 'voxels'}
    "operate_on() -> collections.abc.MutableSet[str]\n\nA set of strings specifying the kinds of data elements this modifier should act on. By default, the set includes all types of data elements supported by the modifier. \n\nDefault: ``{'particles', 'vector_properties', 'cell', 'surfaces', 'dislocations', 'lines', 'vectors', 'voxels'}``"
    reduced_coords: bool = False
    "reduced_coords() -> bool\n\nControls whether the translation vector (fourth column of the :py:attr:`.transformation` matrix) is specified in reduced cell coordinates or in absolute Cartesian coordinates. \n\nIf set to ``False``, the modifier applies the transformation :math:`\\mathbf{x}' =  \\mathbf{M} \\cdot \\mathbf{x} + \\mathbf{t}` to Cartesian input points :math:`\\mathbf{x}`. Here, :math:`\\mathbf{M}` refers to the 3x3 linear part of the affine :py:attr:`.transformation` matrix and :math:`\\mathbf{t}` to the translational part (fourth matrix column). \n\nIf set to ``True``, the modifier applies the transformation :math:`\\mathbf{x}' =  \\mathbf{M} \\cdot (\\mathbf{x} + \\mathbf{H} \\cdot \\mathbf{t})`. Here, :math:`\\mathbf{H}` refers to the 3x3 cell matrix formed by the three edge vectors of the simulation cell. \n\nNote that this option only applies if :py:attr:`relative_mode` is active. \n\nDefault: ``False``"
    relative_mode: bool = True
    'relative_mode() -> bool\n\nSelects the operation mode of the modifier.\n\nIf set to ``True``, the modifier transforms data elements by applying the specified :py:attr:`transformation` matrix.\n\nIf set to ``False``, the modifier determines the effective transformation dynamically based on the current shape of the :py:class:`SimulationCell` and the specified :py:attr:`target_cell` matrix. The old cell will be mapped to the new shape using an appropriate affine transformation. \n\nDefault: ``True``'
    target_cell: ArrayLike = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
    'target_cell() -> numpy.typing.ArrayLike\n\nThis 3x4 matrix specifies the target cell shape. It is only used if :py:attr:`relative_mode` == ``False``. \n\nThe first three columns of the matrix specify the three edge vectors of the target cell. The fourth column specifies the origin vector of the target cell. \n\nThe following code shows how to scale the simulation box, whose shape may vary with simulation time, back to the initial shape at frame 0, including the cell\'s contents. As a result, the output dataset generated by the modifier will have a constant simulation cell size. \n\n```python\n  pipeline = import_file("input/simulation.*.dump")\n  modifier = AffineTransformationModifier(\n          relative_mode = False,\n          target_cell = pipeline.compute(0).cell[...]\n  )\n  pipeline.modifiers.append(modifier)\n```'
    transformation: ArrayLike = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]]
    'transformation() -> numpy.typing.ArrayLike\n\nThe 3x4 transformation matrix to apply to input elements. The first three matrix columns define the linear part of the transformation. The fourth column specifies the translation vector. \n\nNote that this matrix describes a *relative* transformation and is used only if :py:attr:`relative_mode` == ``True``. See the :py:attr:`.reduced_coords` field for a definition of the precise coordinate transformation that is specified by this matrix. \n\nDefault: ``[[1,0,0,0], [0,1,0,0], [0,0,1,0]]``'

@dataclass(kw_only=True)
class AmbientOcclusionModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

Performs a quick lighting calculation to modulate the brightness of particles according to the degree of occlusion by other particles.
See the corresponding user manual page for more information.

This modifier should always be inserted at the end of a data pipeline, after particle filtering and coloring is done. 
It can help to improve the visual appearance of complex particle structures which are rendered with the default
:py:class:`OpenGLRenderer`. It is *not* needed when using the :py:class:`TachyonRenderer`
or :py:class:`OSPRayRenderer`, which both perform the ambient occlusion 
calculation on the fly as part of the rendering process.

  This modifier performs the lighting calculation using the OpenGL interface and GPU accelleration, which requires
  a graphics environment. See the :py:class:`OpenGLRenderer` class for more information
  on how to enable OpenGL graphics for a Python script running in the terminal. 

Inputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Position``
      - The input coordinates of the particles.
    * - ``Color``
      - The original per-particle colors.

Outputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Color``
      - The new per-particle colors, which have been modulated by the occlusion factor computed by the modifier for each particle."""
    buffer_resolution: int = 3
    'buffer_resolution() -> int\n\nA positive integer controlling the resolution of the internal render buffer, which is used to compute how much light each particle receives. For large datasets, where the size of particles is small compared to the simulation dimensions, a higher buffer resolution should be used.\n\n:Valid range: [1, 4]\nDefault: ``3``'
    intensity: float = 0.7
    'intensity() -> float\n\nControls the strength of the shading effect. \n\n:Valid range: [0.0, 1.0]\nDefault: ``0.7``'
    sample_count: int = 40
    'sample_count() -> int\n\nThe number of light exposure samples to compute. More samples give a more even light distribution but take longer to compute.\n\nDefault: ``40``'

    @staticmethod
    def is_available() -> bool:
        ...

@dataclass(kw_only=True)
class AssignColorModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

This modifier assigns a uniform color to all selected data elements.
See also the corresponding user manual page for more information.
The modifier can operate on different data elements:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Data element
      -
    * - ``particles``
      - Sets the ``Color`` property of selected particles.
    * - ``vectors``
      - Sets the ``Vector Color`` property of selected particles.
    * - ``bonds``
      - Sets the ``Color`` property of selected bonds.
    * - ``surface_vertices``
      - Sets the ``Color`` property of selected vertices of a :py:class:`SurfaceMesh` structure.
    * - ``surface_faces``
      - Sets the ``Color`` property of selected faces of a :py:class:`SurfaceMesh` structure.
    * - ``surface_regions``
      - Sets the ``Color`` property of selected spatial regions of a :py:class:`SurfaceMesh` structure.

By default the modifier will act on particles. You can change this by setting the :py:attr:`operate_on` field.
The modifier uses the input :py:class:`Property` named ``Selection`` to determine which subset of data elements should
be assigned the specified :py:attr:`color`. Data elements whose ``Selection`` property is zero, keep their current color.
In case the ``Selection`` property does not exist at all, the modifier assigns the color to every data element.

OVITO provides various modifiers for creating selections: :py:class:`SelectTypeModifier`, :py:class:`ExpressionSelectionModifier`,
:py:class:`SliceModifier`, :py:class:`InvertSelectionModifier`.

Examples

Select some particles with the :py:class:`SelectTypeModifier` and give them a new color:

    ```python
  pipeline.modifiers.append(SelectTypeModifier(types={1,3}))
  pipeline.modifiers.append(AssignColorModifier(color=(1.0, 0.3, 0.3)))
```

Use the :py:class:`CalculateDisplacementsModifier` to calculate the atomic displacement vectors and visualize them
as arrows in rendered images. Particularly large displacement arrows can be highlighted using the :py:class:`AssignColorModifier`:

    ```python
  # Calculate and visualize atomic displacement vectors:
  modifier = CalculateDisplacementsModifier()
  modifier.vis.enabled = True
  pipeline.modifiers.append(modifier)
  
  # Select particles with large displacement magnitudes:
  pipeline.modifiers.append(ExpressionSelectionModifier(expression='DisplacementMagnitude > 1.2'))
  # Highlight large displacement vectors using a special color: 
  pipeline.modifiers.append(AssignColorModifier(operate_on='vectors', color=(0.0, 0.9, 1.0)))
```"""
    color: ovito.vis.Color = (0.3, 0.3, 1.0)
    'color() -> tuple[float, float, float]\n\nThe RGB color that will be assigned to all selected elements by the modifier.\n\nDefault: ``(0.3, 0.3, 1.0)``'
    operate_on: str = 'particles'
    "operate_on() -> str\n\nSelects the kind of data elements this modifier should operate on. Supported values are: ``'particles'``, ``'bonds'``, ``'vectors'``, ``'surface_vertices'``, ``'surface_faces'``, ``'surface_regions'``. \n\nDefault: ``'particles'``"

@dataclass(kw_only=True)
class AtomicStrainModifier(ovito.pipeline.ReferenceConfigurationModifier):
    """Base: :py:class:`ovito.pipeline.ReferenceConfigurationModifier`

Computes the atomic-level deformation with respect to a reference configuration.
See the corresponding user manual page
for more information.

The modifier is a subclass of :py:class:`ReferenceConfigurationModifier`, which provides
the programming interface for specifying the reference configuration and how particle displacements get calculated.
By default, frame 0 of the processed simulation sequence is used as static reference configuration.

Outputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Shear Strain``
      - The von Mises shear strain invariant of the computed atomic Green-Lagrangian strain tensor.
    * - ``Volumetric Strain``
      - One third of the trace of the computed atomic Green-Lagrangian strain tensor.
    * - ``Strain Tensor``
      - The six components of the symmetric Green-Lagrangian strain tensor. Only if :py:attr:`output_strain_tensors` was set to ``True``.
    * - ``Deformation Gradient``
      - The nine components of the atomic deformation gradient tensor. Only if :py:attr:`output_deformation_gradients` was set to ``True``.
    * - ``Stretch Tensor``
      - The six components of the symmetric right stretch tensor U in the polar decomposition F=RU. Only if :py:attr:`output_stretch_tensors` was set to ``True``.
    * - ``Rotation``
      - The atomic microrotation obtained from the polar decomposition F=RU as a 4-component quaternion. Only if :py:attr:`output_rotations` was set to ``True``.
    * - ``Nonaffine Squared Displacement``
      - The D\\ :sup:`2`\\ :sub:`min` measure of Falk & Langer, which describes the non-affine part of the local deformation. Only if :py:attr:`output_nonaffine_squared_displacements` was set to ``True``.
    * - ``Selection``
      -  Set to a non-zero value for particles for which the modifier failed to determine a local deformation tensor, because they do not have enough neighbors
         within the specified :py:attr:`cutoff` distance. Only if :py:attr:`select_invalid_particles` was set to ``True``.
         The selected particles without valid deformation values can subsequently be removed using a :py:class:`DeleteSelectedModifier`.

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Global attributes
      -
    * - ``AtomicStrain.invalid_particle_count``
      - Number of particles for which the modifier could not compute a deformation tensor, because they do not have enough neighbors
        within the specified :py:attr:`cutoff` distance. You typically should increase the cutoff distance if this value is non-zero."""
    cutoff: float = 3.0
    'cutoff() -> float\n\nThe spatial range up to which neighboring atoms will be taken into account to calculate the local strain measure.\n\nDefault: ``3.0``'
    output_deformation_gradients: bool = False
    'output_deformation_gradients() -> bool\n\nControls the output of the per-particle deformation gradient tensors. If ``False``, the computed tensors are not output as a particle property to save memory.\n\nDefault: ``False``'
    output_nonaffine_squared_displacements: bool = False
    'output_nonaffine_squared_displacements() -> bool\n\nEnables the computation of the squared magnitude of the non-affine part of the atomic displacements. The computed values are output in the ``"Nonaffine Squared Displacement"`` particle property.\n\nDefault: ``False``'
    output_rotations: bool = False
    'output_rotations() -> bool\n\nControls the calculation of the per-particle rotations.\n\nDefault: ``False``'
    output_strain_tensors: bool = False
    'output_strain_tensors() -> bool\n\nControls the output of the per-particle strain tensors. If ``False``, the computed strain tensors are not output as a particle property to save memory.\n\nDefault: ``False``'
    output_stretch_tensors: bool = False
    'output_stretch_tensors() -> bool\n\nControls the calculation of the per-particle stretch tensors.\n\nDefault: ``False``'
    select_invalid_particles: bool = True
    'select_invalid_particles() -> bool\n\nIf ``True``, the modifier selects the particle for which the local strain tensor could not be computed (because of an insufficient number of neighbors within the cutoff).\n\nDefault: ``True``'

@dataclass(kw_only=True)
class BondAnalysisModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

Computes the bond angle distribution and the bond length distribution for the particle system with bonds.
See the corresponding user manual page for more information on this modifier.

The input for this modifier must contain a :py:class:`Bonds` data object with the bond topology information,
which may either be loaded from the simulation file or by first inserting the :py:class:`CreateBondsModifier`
into the pipeline.

The modifier outputs the computed bond length and bond angle histograms as two :py:class:`DataTable` objects,
which may be exported to an output text file using the :py:func:`export_file` function or retrieved
from the pipeline's output data collection, see the `DataCollection.tables` dictionary
and the code examples following below.

Inputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Position``
      - The coordinates of the input particles.
    * - ``Particle Type``
      - Required if :py:attr:`partition` is set to ``ByParticleType``.
    * - ``Selection``
      - Required if :py:attr:`partition` is set to ``ByParticleSelection``.

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Bond properties
      -
    * - ``Topology``
      - The input bonds list.
    * - ``Bond Type``
      - Required if :py:attr:`partition` is set to ``ByBondType``.
    * - ``Selection``
      - Required if :py:attr:`partition` is set to ``ByBondSelection``.

Outputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - :py:class:`DataTable`
      -
    * - ``bond-angle-distr``
      - The bond angle distribution histogram computed by the modifier.
    * - ``bond-length-distr``
      - The bond length distribution histogram computed by the modifier.

Examples:

The following Python script demonstrates how to load a particle system (without bond topology) from a file,
create the bond topology using OVITO's :py:class:`CreateBondsModifier`, and then compute
the bond angle distribution. Finally, the histogram is written to an output text file using the :py:func:`export_file` function.

```python
  from ovito.io import import_file, export_file
  from ovito.modifiers import BondAnalysisModifier, CreateBondsModifier
  
  # Set up data pipeline:
  pipeline = import_file('input/simulation.dump')
  pipeline.modifiers.append(CreateBondsModifier(cutoff = 3.2))
  pipeline.modifiers.append(BondAnalysisModifier(bins = 100))
  
  # Export bond angle distribution to an output text file.
  export_file(pipeline, 'output/bond_angles.txt', 'txt/table', key='bond-angle-distr')
  
  # Convert bond length histogram to a NumPy array and print it to the terminal.
  data = pipeline.compute()
  print(data.tables['bond-length-distr'].xy())
```

The script above computes the instantaneous distribution for the initial frame in the simulation file only.
To compute an average bond angle distribution for the entire MD trajectory, you can make use of the
:py:class:`TimeAveragingModifier`:

```python
  from ovito.io import *
  from ovito.modifiers import *
  
  # Load a LAMMPS topology data file containing bond information:
  pipeline = import_file('input/input.data', atom_style='bond')
  
  # Load trajectory data from separate LAMMPS dump file.
  traj = LoadTrajectoryModifier()
  traj.source.load('input/trajectory.dump')
  pipeline.modifiers.append(traj)
  
  # Calculate instantaneous bond angle distribution. 
  pipeline.modifiers.append(BondAnalysisModifier())
  
  # Perform time averaging of the DataTable 'bond-angle-distr'.
  pipeline.modifiers.append(TimeAveragingModifier(operate_on='table:bond-angle-distr'))
  
  # Compute and export the time-averaged histogram to a text file.
  export_file(pipeline, "output/avg_bond_angles.txt", 'txt/table', key='bond-angle-distr[average]')
```"""

    class Partition(enum.Enum):
        """"""
        Off = enum.auto()
        ByBondType = enum.auto()
        ByBondSelection = enum.auto()
        ByParticleType = enum.auto()
        ByParticleSelection = enum.auto()
    partition: BondAnalysisModifier.Partition = Partition.Off
    'partition() -> BondAnalysisModifier.Partition\n\nThis mode parameter instructs the modifier to split the bond length and bond angle distributions into partial distributions, one for each particle or bond type combination. Available partitioning modes are: \n\n   * ``BondAnalysisModifier.Partition.Off``\n   * ``BondAnalysisModifier.Partition.ByBondType``\n   * ``BondAnalysisModifier.Partition.ByBondSelection``\n   * ``BondAnalysisModifier.Partition.ByParticleType``\n   * ``BondAnalysisModifier.Partition.ByParticleSelection``\n\n\nIf partitioning is turned on, the :py:attr:`y` property array of the two :py:class:`DataTable` histograms computed by the modifier will have multiple vector components, one for each partial distribution. The follow code example shows how to access the partial histograms: \n\n```python\n  from ovito.io import import_file\n  from ovito.modifiers import BondAnalysisModifier\n  \n  # Load LAMMPS data file with bonds and compute bond angle distribution partitioned by bond type:\n  pipeline = import_file(\'input/input.data\', atom_style=\'bond\')\n  pipeline.modifiers.append(BondAnalysisModifier(partition=BondAnalysisModifier.Partition.ByBondType))\n  data = pipeline.compute()\n  \n  # Retrieve y-property array from output DataTable:\n  histogram = data.tables[\'bond-angle-distr\'].y\n  \n  # Print the individual columns of the vector-valued property array:\n  for column, name in enumerate(histogram.component_names):\n      print("Angle distribution for bond types:", name)\n      print(histogram[:,column])\n```\n\nDefault: ``BondAnalysisModifier.Partition.Off``'
    bins: int = 180
    'bins() -> int\n\nThe number of bins in the length and angle histograms generated by the modifier. \n\nDefault: ``180``'
    cosine_mode: bool = False
    'cosine_mode() -> bool\n\nIf set to ``True``, the modifier will calculate the distribution of the cosines of the bond angles instead of the angles themselves, and the bond angle histogram output by the modifier will then extend over the value interval [-1, +1] instead of [0, 180] degrees. \n\nDefault: ``False``'
    length_cutoff: float = 4.0
    'length_cutoff() -> float\n\nMaximum bond length at which the bond length distribution gets truncated. Together with the :py:attr:`bins` parameter, this value determines the bin size of the bond length histogram output by the modifier. \n\nNote that this parameter does not affect the computation of the bond angle distribution. All bonds are included in the angle distribution, even if their length exceeds the cutoff. \n\nDefault: ``4.0``'

@dataclass(kw_only=True)
class CalculateDisplacementsModifier(ovito.pipeline.ReferenceConfigurationModifier):
    """Base: :py:class:`ovito.pipeline.ReferenceConfigurationModifier`

Computes the displacement vectors of particles with respect to a reference configuration.
See the corresponding user manual page
for more information.

The modifier is a subclass of :py:class:`ReferenceConfigurationModifier`, which provides
the programming interface for specifying the reference configuration and how particle displacements get calculated.
By default, frame 0 of the processed simulation sequence is used as static reference configuration.

Outputs:

.. list-table::
    :widths: 35 65
    :header-rows: 1

    * - Particle properties
      -
    * - ``Displacement``
      - The computed displacement vectors.
    * - ``Displacement Magnitude``
      - The length of the computed displacement vectors."""
    vis: ovito.vis.VectorVis = ovito.vis.VectorVis(enabled=False, alignment=ovito.vis.VectorVis.Alignment.Head, title='Displacements')
    'vis() -> ovito.vis.VectorVis\n\nA :py:class:`VectorVis` element controlling the visual representation of the computed displacement vectors. Note that the computed displacement vectors are hidden by default. You can enable the visualization of arrows as follows: \n\n```python\n  modifier = CalculateDisplacementsModifier()\n  modifier.vis.enabled = True\n  modifier.vis.color = (0.8, 0.0, 0.5)\n```'

@dataclass(kw_only=True)
class CentroSymmetryModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

Computes the centrosymmetry parameter (CSP) of each particle, which is a measure of the local lattice disorder
around a particle in centrosymmetric crystal lattices.
See the corresponding user manual page
for more information.

Inputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Position``
      - The input coordinates of the particles.

Outputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Centrosymmetry``
      - The non-negative CSP value computed for each particle. Values close to zero mean neighboring
        particles are in a perfect centrosymmetric arrangement."""

    class Mode(enum.Enum):
        """"""
        Conventional = enum.auto()
        Matching = enum.auto()
    mode: CentroSymmetryModifier.Mode = Mode.Conventional
    'mode() -> CentroSymmetryModifier.Mode\n\nSelects how pair-wise opposite neighbors should be picked during the CSP calculation. Valid modes are:\n\n  * ``CentroSymmetryModifier.Mode.Conventional``\n  * ``CentroSymmetryModifier.Mode.Matching``\n\n\nSee the user manual for more information on these two modes. \n\nDefault: ``CentroSymmetryModifier.Mode.Conventional``'
    num_neighbors: int = 12
    'num_neighbors() -> int\n\nThe number of nearest neighbors to take into account for the computation, e.g. \n\n  * 12 for FCC crystals\n  * 8 for BCC crystals\n\n\nDefault: ``12``'
    only_selected: bool = False
    "only_selected() -> bool\n\nSet this to ``True`` to perform the analysis only on the sub-set of currently selected particles. Particles that are *not* selected will be treated as if they did not exist. That means they won't be considered in the centrosymmetry calculation of the surrounding particles and their own centrosymmetry value will be zero. Use a :py:class:`SelectTypeModifier` in your pipeline, for example, to restrict the centrosymmetry analysis to a sub-lattice formed by one species of particles in a multi-component system. \n\nDefault: ``False``"

@dataclass(kw_only=True)
class ChillPlusModifier(StructureIdentificationModifier):
    """Base: :py:class:`ovito.modifiers.StructureIdentificationModifier`

Analyzes the local neighborhood of each particle to identify different structural arrangements of water molecules.
The structure identification is performed using the CHILL+ algorithm.
See the corresponding user manual page
for more information.

Note that this class inherits several important parameter fields from its :py:class:`StructureIdentificationModifier`
base class.

Modifier inputs:

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Particle properties
      -
    * - ``Position``
      - The input coordinates of the particles.
    * - ``Selection``
      - The selection state of the input particles. Only needed if :py:attr:`~StructureIdentificationModifier.only_selected` is ``True``.

Modifier outputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Structure Type``
      - The structure type computed by the algorithm for each particle, encoded as an integer value:

        ============= =========================================================
        Numeric id    Python constant
        ============= =========================================================
        0             ``ChillPlusModifier.Type.OTHER``
        1             ``ChillPlusModifier.Type.HEXAGONAL_ICE``
        2             ``ChillPlusModifier.Type.CUBIC_ICE``
        3             ``ChillPlusModifier.Type.INTERFACIAL_ICE``
        4             ``ChillPlusModifier.Type.HYDRATE``
        5             ``ChillPlusModifier.Type.INTERFACIAL_HYDRATE``
        ============= =========================================================
    * - ``Color``
      - Particle coloring to indicate the identified structure type for each particle; only if :py:attr:`~StructureIdentificationModifier.color_by_type` is ``True``.
        See the :py:attr:`~StructureIdentificationModifier.structures` array on how to customize the colors.

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Global attributes
      -
    * - ``ChillPlus.counts.OTHER``
      - Number of particles not matching any of the known structure types.
    * - ``ChillPlus.counts.HEXAGONAL_ICE``
      - Number of particles identified as hexagonal ice.
    * - ``ChillPlus.counts.CUBIC_ICE``
      - Number of particles identified as cubic ice.
    * - ``ChillPlus.counts.INTERFACIAL_ICE``
      - Number of particles identified as interfacial ice.
    * - ``ChillPlus.counts.HYDRATE``
      - Number of particles identified as hydrate.
    * - ``ChillPlus.counts.INTERFACIAL_HYDRATE``
      - Number of particles identified as interfacial hydrate.

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Data tables
      -
    * - ``structures``
      - A bar chart with the particle counts for each structure type identified by the modifier.
        You can retrieve this :py:class:`DataTable` from the `DataCollection.tables` dictionary."""

    class Type(enum.IntEnum):
        """"""
        OTHER = enum.auto()
        HEXAGONAL_ICE = enum.auto()
        CUBIC_ICE = enum.auto()
        INTERFACIAL_ICE = enum.auto()
        HYDRATE = enum.auto()
        INTERFACIAL_HYDRATE = enum.auto()
    cutoff: float = 3.5
    'cutoff() -> float\n\nThe cutoff distance for bonds between water molecules. \n\nDefault: ``3.5``'

@dataclass(kw_only=True)
class ClearSelectionModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

This modifier clears the current selection by removing the ``Selection`` property from a :py:class:`PropertyContainer`
such that subsequent modifiers in the pipeline won't see it.
See also the corresponding user manual page for more information.
The modifier can operate on different kinds of data elements:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Data element
      -
    * - ``particles``
      - Removes the ``Selection`` property of particles.
    * - ``bonds``
      - Removes the ``Selection`` property of bonds.
    * - ``voxels``
      - Removes the ``Selection`` property of voxel grid cells.

By default the modifier will act on particles. You can change this by setting the :py:attr:`operate_on` field."""
    operate_on: str = 'particles'
    "operate_on() -> str\n\nSelects the kind of data elements this modifier should operate on. Supported values are: ``'particles'``, ``'bonds'``, ``'voxels'``. \n\nDefault: ``'particles'``"

@dataclass(kw_only=True)
class ClusterAnalysisModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

This modifier groups particles into disconnected clusters based on the chosen connectivity criterion.
See the corresponding user manual page
for more information.

Outputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Cluster``
      - The 1-based numeric ID of the cluster each particle was assigned to by the modifier.

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Global attributes
      -
    * - ``ClusterAnalysis.cluster_count``
      - Total number of clusters produced by the modifier. Cluster IDs range from 1 to this number.
    * - ``ClusterAnalysis.largest_size``
      - Number of particles in the largest cluster (cluster ID 1). Only computed if :py:attr:`sort_by_size` is set to ``True``.

.. list-table::
    :widths: 15 85
    :header-rows: 1

    * - :py:class:`DataTable`
      -
    * - ``clusters``
      - The list of clusters identified by the modifier.
        You can retrieve this :py:class:`DataTable` from the `DataCollection.tables` dictionary.
        It consists of several data columns, each being a separate :py:class:`Property` object:

          * ``Cluster Identifier``: The numeric ID of each identified cluster (1-based).
          * ``Cluster Size``: The number of particles in each identified cluster.
          * ``Center of Mass``: The XYZ coordinates of the center of mass of each cluster (only if :py:attr:`compute_com` is set).
          * ``Radius of Gyration``: The radius of gyration of each cluster (only if :py:attr:`compute_gyration` is set) in simulation units of length.
          * ``Gyration Tensor``: The gyration tensor of each cluster (only if :py:attr:`compute_gyration` is set) in simulation units of length squared.
            The tensors are stored as vectors with six components [XX, YY, ZZ, XY, XZ, YZ].

Example:

The following example code demonstrates how to export the data table generated by the modifier to a text file and how to access the
information from Python.

```python
  from ovito.io import import_file, export_file
  from ovito.modifiers import ClusterAnalysisModifier
  import numpy
  
  pipeline = import_file("input/simulation.dump")
  pipeline.modifiers.append(ClusterAnalysisModifier(
      cutoff=2.8, 
      sort_by_size=True, 
      compute_com=True))
  
  # Export results of the clustering algorithm to a text file:
  export_file(pipeline, 'output/clusters.txt', 'txt/table', key='clusters')
  
  # Directly access information stored in the DataTable:
  data = pipeline.compute()
  cluster_table = data.tables['clusters']
  print("Centers of mass:")
  print(cluster_table['Center of Mass'][...])
```"""

    class NeighborMode(enum.Enum):
        """"""
        CutoffRange = enum.auto()
        Bonding = enum.auto()
    cluster_coloring: bool = False
    'cluster_coloring() -> bool\n\nEnables the coloring of particles based on their assignment to clusters. Each cluster is represented by a unique random color. \n\nDefault: ``False``'
    compute_com: bool = False
    'compute_com() -> bool\n\nEnables the computation of the center of mass of each cluster. The center coordinates will be output as an extra column named ``Center of Mass`` in the ``clusters`` data table. \n\nDefault: ``False``'
    compute_gyration: bool = False
    'compute_gyration() -> bool\n\nEnables the computation of the radius of gyration and the gyration tensor of each cluster. Both quantities will be output as auxiliary properties to the ``clusters`` data table, see above. \n\nDefault: ``False``'
    cutoff: float = 3.2
    'cutoff() -> float\n\nThe cutoff distance used by the algorithm to form clusters of connected particles. This parameter is only used when :py:attr:`neighbor_mode` is set to ``CutoffRange``; otherwise it is ignored. \n\nDefault: ``3.2``'
    neighbor_mode: ClusterAnalysisModifier.NeighborMode = NeighborMode.CutoffRange
    'neighbor_mode() -> ClusterAnalysisModifier.NeighborMode\n\nSelects the neighboring criterion for the clustering algorithm. Valid values are: \n\n  * ``ClusterAnalysisModifier.NeighborMode.CutoffRange``\n  * ``ClusterAnalysisModifier.NeighborMode.Bonding``\n\n\nIn the first mode (``CutoffRange``), the clustering algorithm treats pairs of particles as neighbors which are within a certain range of each other given by the parameter :py:attr:`cutoff`. \n\nIn the second mode (``Bonding``), particles which are connected by bonds are combined into clusters. Bonds between particles can either be loaded from the input simulation file or dynamically created using for example the :py:class:`CreateBondsModifier` or the :py:class:`VoronoiAnalysisModifier`. \n\nDefault: ``ClusterAnalysisModifier.NeighborMode.CutoffRange``'
    only_selected: bool = False
    'only_selected() -> bool\n\nLets the modifier perform the analysis only for selected particles. In this case, particles which are not selected are treated as if they did not exist and will be assigned cluster ID 0. \n\nDefault: ``False``'
    sort_by_size: bool = False
    'sort_by_size() -> bool\n\nEnables the sorting of clusters by size (in descending order). Cluster 1 will be the largest cluster, cluster 2 the second largest, and so on.\n\nDefault: ``False``'
    unwrap_particles: bool = False
    'unwrap_particles() -> bool\n\nEnables the "unwrapping" of particle coordinates in order to make clusters contiguous that are cut by a periodic simulation cell boundary. \n\nDefault: ``False``'

@dataclass(kw_only=True)
class ColorByTypeModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

Assigns colors to data elements to visualize the discrete per-element values of a typed property, e.g. the residue type or structural type of particles.
See also the corresponding user manual page for more information.
The modifier can operate on different data elements:

.. list-table::
    :widths: 20 80
    :header-rows: 1

    * - Data elements
      -
    * - ``particles``
      - Colors :py:class:`Particles` according to a typed particle property.
    * - ``bonds``
      - Colors :py:class:`Bonds` according to a typed bond property.
    * - ``voxels``
      - Colors :py:class:`VoxelGrid` cells according to a typed property.

By default the modifier will act on particles. You can change this by setting the :py:attr:`operate_on` field.

The modifier's :py:attr:`property` field names the input property according to which the data elements will be colored.
It must be a so-called *typed property*, which means the :py:class:`Property` object has a list of :py:class:`ElementType`
instances attached in its :py:attr:`types` field, defining the mapping of integer values in the property array to element type colors.
Each :py:class:`ElementType` instance associates one numeric :py:attr:`id` with
a corresponding :py:attr:`color`.

Output:

The modifier sets the ``Color`` property of the data elements it operates on. When rendering a picture,
this property determines the visual color of the individual objects.

Examples:

```python
  from ovito.io import import_file
  from ovito.modifiers import ColorByTypeModifier
  
  # Load a GROMACS file, which contains the 'Residue Type' particle property.
  pipeline = import_file("input/1AKI.gro")
  pipeline.add_to_scene()
  
  # Print a table of all residue types defined in the dataset:
  for residue_type in pipeline.compute().particles['Residue Type'].types:
      print(residue_type.id, residue_type.name, residue_type.color)
  
  # Apply the ColorByTypeModifier, giving particles a color based on 
  # the value of their 'Residue Type' property:
  pipeline.modifiers.append(ColorByTypeModifier(
      operate_on='particles', 
      property='Residue Type'))
```

One way to change the colors of the element types is to prepend the :py:class:`ColorByTypeModifier` in the pipeline with a user-defined modifier function
configuring the :py:class:`ElementType` instances that are attached to the typed :py:class:`Property`:

```python
  def setup_residue_colors(frame, data):
      residues = data.particles_['Residue Type_']
      residues.type_by_name_('LYS').color = (1.0, 0.0, 1.0)
      residues.type_by_name_('GLU').color = (0.0, 0.5, 1.0)
  pipeline.modifiers.insert(0, setup_residue_colors)
```

Implementation:

The :py:class:`ColorByTypeModifier` is functionally equivalent to (but more efficient than) the following user-defined modifier function:

```python
  def color_by_type(frame, data, property_name='Residue Type'):
      input_property = data.particles[property_name]
      output_colors = data.particles_.create_property('Color')
      for index, type_id in enumerate(input_property):
          element_type = input_property.type_by_id(type_id)
          output_colors[index] = element_type.color
```"""
    clear_selection: bool = True
    'clear_selection() -> bool\n\nControls whether the current element selection is cleared to reveal the assigned colors in the interactive viewports of OVITO, which may otherwise be masked by the highlighting of selected elements. This option only has an effect if the :py:attr:`.only_selected` option is also turned on. Otherwise the modifier never clears the current element selection. \n\nDefault: ``True``'
    only_selected: bool = False
    'only_selected() -> bool\n\nIf ``True``, only selected data elements will be assigned a color by the modifier and the colors of unselected elements will be preserved; if ``False``, all data elements will be colored.\n\nDefault: ``False``'
    operate_on: str = 'particles'
    "operate_on() -> str\n\nControls the kind of data elements that are being colored by this modifier. Supported values are: ``'particles'``, ``'bonds'``, ``'voxels'``. \n\nDefault: ``'particles'``"
    property: str = 'Particle Type'
    "property() -> str\n\nThe name of the typed property to use as input; must be an integer property with attached :py:class:`ElementType` instances. \n\nFor particles, typical input properties are ``'Particle Type'``, ``'Residue Type'`` or ``'Structure Type'``. When coloring bonds, the ``'Bond Type'`` property is an example for a typed property. \n\nDefault: ``'Particle Type'``"

@dataclass(kw_only=True)
class ColorCodingModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

Colors elements to visualize one of their properties.
See also the corresponding user manual page for more information.
The modifier can operate on different data elements:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Data element
      -
    * - ``particles``
      - Assigns a local color to particles by setting their ``Color`` property.
    * - ``vectors``
      - Colors the vector arrows associated with particles by setting their ``Vector Color`` property.
    * - ``bonds``
      - Assigns a local color to bonds by setting their ``Color`` property.
    * - ``voxels``
      - Colors the cells of a :py:class:`VoxelGrid`.
    * - ``surface_vertices``
      - Colors the vertices of a :py:class:`SurfaceMesh` object.
    * - ``surface_faces``
      - Colors the faces of a :py:class:`SurfaceMesh` object.
    * - ``surface_regions``
      - Colors the volumetric regions of a :py:class:`SurfaceMesh` object.
    * - ``lines``
      - Assigns colors to the vertices of a :py:class:`Lines` object to visualize a local property.

By default the modifier will act on particles. You can change this by setting the :py:attr:`operate_on` field.

Example:

```python
  pipeline.modifiers.append(ColorCodingModifier(
      property = 'Potential Energy',
      gradient = ColorCodingModifier.Hot()
  ))
```

If, as in the example above, the :py:attr:`start_value` and :py:attr:`end_value` parameters are not explicitly specified,
then modifier will automatically adjust the mapping interval to fully cover the range of input property values
(dynamically for each trajectory frame).

The :py:class:`ColorLegendOverlay` may be used in conjunction with a :py:class:`ColorCodingModifier`
to include a color legend in rendered images."""

    class Gradient:

        def __init__(self, color_array: ArrayLike) -> None:
            """Initialize self.  See help(type(self)) for accurate signature."""
            ...

    class BlueWhiteRed(Gradient):

        def __init__(self) -> None:
            ...

    class CyclicRainbow(Gradient):

        def __init__(self) -> None:
            ...

    class Grayscale(Gradient):

        def __init__(self) -> None:
            ...

    class Hot(Gradient):

        def __init__(self) -> None:
            ...

    class Jet(Gradient):

        def __init__(self) -> None:
            ...

    class Magma(Gradient):

        def __init__(self) -> None:
            ...

    class Rainbow(Gradient):

        def __init__(self) -> None:
            ...

    class Viridis(Gradient):

        def __init__(self) -> None:
            ...

    class Image(Gradient):

        def __init__(self, imagefile: str) -> None:
            ...
    auto_adjust_range: bool = True
    "auto_adjust_range() -> bool\n\nControls the automatic adjustment of the modifier's mapping interval to always cover the full range of input values. If turned on (the default), the modifier adaptively adjusts the value-to-color mapping to the current min/max range of input values (at the current timestep). If turned off, the mapping is fixed and is determined by the :py:attr:`start_value` and :py:attr:`end_value` parameters of the modifier. \n\nSetting the :py:attr:`start_value` or :py:attr:`end_value` parameters to some value implicitly changes :py:attr:`auto_adjust_range` to ``False``. Furthermore, note that the automatic range is always determined over the complete set of input elements irrespective of the :py:attr:`only_selected` option. \n\nDefault: ``True``"
    symmetric_range: bool = False
    'symmetric_range() -> bool\n\nWhen this option is enabled, the modifier will automatically adjust the value range to be symmetric around 0. Manually setting the maximum value will also adjust the minimum value to match.\n\nDefault: ``False``'
    end_value: float = 0.0
    'end_value() -> float\n\nThis parameter determines, together with :py:attr:`start_value`, the linear mapping of input property values to colors. It is only used if :py:attr:`auto_adjust_range` is turned off, which happens automatically as soon as you assign some value to this modifier parameter. \n\nDefault: ``0.0``'
    gradient: ColorCodingModifier.Gradient = Rainbow()
    'The color gradient used to map normalized property values to colors. Available gradient types are:\n\n * ``ColorCodingModifier.BlueWhiteRed()``\n * ``ColorCodingModifier.CyclicRainbow()``\n * ``ColorCodingModifier.Grayscale()``\n * ``ColorCodingModifier.Hot()``\n * ``ColorCodingModifier.Jet()``\n * ``ColorCodingModifier.Magma()``\n * ``ColorCodingModifier.Rainbow()`` (default)\n * ``ColorCodingModifier.Viridis()``\n * ``ColorCodingModifier.Gradient(<(N,3) array>)``\n * ``ColorCodingModifier.Image("<image file path>")``\n\nThe ``Gradient`` constructor lets you define your own coloring scheme and takes an array of dimensions *N* x 3 containing a table of colors (RGB values in the range [0-1]). The color coding modifier will linearly interpolate between the *N* colors of the table. \nThe ``Image`` constructor expects the path to an image file on disk, which will be used to create a custom color gradient from the first row of pixels in the image.\n\nCode example:\n\n```python\n  color_table = [\n      (1.0, 0.0, 0.0),  # red\n      (1.0, 1.0, 0.0),  # yellow\n      (1.0, 1.0, 1.0)   # white\n  ]\n  modifier = ColorCodingModifier(\n      property = \'Position.X\',\n      gradient = ColorCodingModifier.Gradient(color_table)\n  )\n```'
    only_selected: bool = False
    'only_selected() -> bool\n\nIf ``True``, only selected elements will be affected by the modifier and the existing colors of unselected elements will be preserved; if ``False``, all elements will be colored.\n\nDefault: ``False``'
    operate_on: str = 'particles'
    "operate_on() -> str\n\nSelects the kind of data elements this modifier should operate on. Supported values are: ``'particles'``, ``'bonds'``, ``'vectors'``, ``'voxels'``, ``'surface_vertices'``, ``'surface_faces'``, ``'surface_regions'``, ``'lines'``. \n\nDefault: ``'particles'``"
    property: str = ''
    'property() -> str\n\nThe name of the input property that should be used to color elements. \n\nWhen the input property has multiple components, then a component name must be appended to the property base name, e.g. ``"Velocity.X"``.'
    start_value: float = 0.0
    'start_value() -> float\n\nThis parameter determines, together with :py:attr:`end_value`, the linear mapping of input property values to colors. It is only used if :py:attr:`auto_adjust_range` is turned off, which happens automatically as soon as you assign some value to this modifier parameter. \n\nDefault: ``0.0``'

@dataclass(kw_only=True)
class CombineDatasetsModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

This modifier loads a set of particles from a separate simulation file and merges them into the primary dataset.
See also the corresponding user manual page for more information.

Example:

```python
  from ovito.io import import_file, export_file
  from ovito.modifiers import CombineDatasetsModifier
  
  # Load a first set of particles.
  pipeline = import_file('input/first_file.dump')
  
  # Insert the particles from a second file into the dataset. 
  modifier = CombineDatasetsModifier()
  modifier.source.load('input/second_file.dump')
  pipeline.modifiers.append(modifier)
  
  # Export combined dataset to a new file.
  export_file(pipeline, 'output/combined.dump', 'lammps/dump',
              columns = ['Position.X', 'Position.Y', 'Position.Z'])
```"""
    source: Optional[ovito.pipeline.FileSource] = ovito.pipeline.FileSource()
    'source() -> ovito.pipeline.FileSource\n\n\nA :py:class:`FileSource` that provides the set of particles to be merged. You can call its :py:meth:`load` function to load a data file as shown in the code example above.'

@dataclass(kw_only=True)
class CommonNeighborAnalysisModifier(StructureIdentificationModifier):
    """Base: :py:class:`ovito.modifiers.StructureIdentificationModifier`

This modifier analyzes the local neighborhood of each particle to identify simple crystalline structures.
The structure identification is performed using the Common Neighbor Analysis (CNA) method.
See the corresponding user manual page
for more information on this modifier and the structure identification algorithm it implements.

Note that this class inherits several important parameter fields from its :py:class:`StructureIdentificationModifier`
base class.

Modifier inputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Position``
      - The input coordinates of the particles. They are not used if :py:attr:`mode` is set to ``BondBased``.
    * - ``Selection``
      - The selection state of the input particles. Only needed if :py:attr:`~StructureIdentificationModifier.only_selected` is set to ``True``.

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Bond properties
      -
    * - ``Topology``
      - The input bond topology; required only if :py:attr:`mode` is set to ``BondBased``.

Modifier outputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Structure Type``
      - The structure type computed by the algorithm for each particle, encoded as an integer value:

        ============= =========================================================
        Numeric id    Python constant
        ============= =========================================================
        0             ``CommonNeighborAnalysisModifier.Type.OTHER``
        1             ``CommonNeighborAnalysisModifier.Type.FCC``
        2             ``CommonNeighborAnalysisModifier.Type.HCP``
        3             ``CommonNeighborAnalysisModifier.Type.BCC``
        4             ``CommonNeighborAnalysisModifier.Type.ICO``
        ============= =========================================================
    * - ``Color``
      - Particle coloring to indicate the identified structure type for each particle; only if :py:attr:`~StructureIdentificationModifier.color_by_type` is ``True``.
        See the :py:attr:`~StructureIdentificationModifier.structures` array on how to customize the colors.

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Global attributes
      -
    * - ``CommonNeighborAnalysis.counts.OTHER``
      - Number of particles not matching any of the recognized structure types.
    * - ``CommonNeighborAnalysis.counts.FCC``
      - Number of particles identified as face-centered cubic.
    * - ``CommonNeighborAnalysis.counts.HCP``
      - Number of particles identified as hexagonal close packed.
    * - ``CommonNeighborAnalysis.counts.BCC``
      - Number of particles identified as body-centered cubic.
    * - ``CommonNeighborAnalysis.counts.ICO``
      - Number of particles identified as icosahedral.

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Data tables
      -
    * - ``structures``
      - A bar chart with the particle counts for each structure type identified by the modifier.
        You can retrieve this :py:class:`DataTable` from the `DataCollection.tables` dictionary."""

    class Type(enum.IntEnum):
        """"""
        OTHER = enum.auto()
        FCC = enum.auto()
        HCP = enum.auto()
        BCC = enum.auto()
        ICO = enum.auto()

    class Mode(enum.Enum):
        """"""
        FixedCutoff = enum.auto()
        AdaptiveCutoff = enum.auto()
        IntervalCutoff = enum.auto()
        BondBased = enum.auto()
    cutoff: float = 3.2
    'cutoff() -> float\n\nThe cutoff radius used for the conventional common neighbor analysis. This parameter is only used if :py:attr:`mode` is ``CommonNeighborAnalysisModifier.Mode.FixedCutoff``.\n\nDefault: ``3.2``'
    mode: CommonNeighborAnalysisModifier.Mode = Mode.AdaptiveCutoff
    'mode() -> CommonNeighborAnalysisModifier.Mode\n\nSelects the algorithm type. One of the following constants:\n\n  * ``CommonNeighborAnalysisModifier.Mode.FixedCutoff``\n  * ``CommonNeighborAnalysisModifier.Mode.AdaptiveCutoff``\n  * ``CommonNeighborAnalysisModifier.Mode.IntervalCutoff``\n  * ``CommonNeighborAnalysisModifier.Mode.BondBased``\n\n\nDefault: ``CommonNeighborAnalysisModifier.Mode.AdaptiveCutoff``'

@dataclass(kw_only=True)
class ComputePropertyModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

Evaluates a user-defined math expression for every input element and stores the results in an output property.
See also the corresponding user manual page for more information.

The modifier can operate on different kinds of data elements, such as particles, bonds, or voxels,
which is controlled by the :py:attr:`operate_on` field.

Example:

```python
  from ovito.modifiers import ComputePropertyModifier
  
  pipeline.modifiers.append(ComputePropertyModifier(
      output_property = 'Color',
      expressions = ['Position.X / CellSize.X', '0.0', '0.5']
  ))
```

In many cases you can replace the use of a :py:class:`ComputePropertyModifier` in a pipeline with an equivalent Python modifier function
  performing the same task. The latter may be the better choice if the computation involves complex indexing operations or conditions,
  which can be accomplished more easily in the Python language or using NumPy array computations.

  For example, the :py:class:`ComputePropertyModifier` creating the ``Color`` particle property in the code snipped above can be
  replaced by the following equivalent user-defined modifier function:

    ```python
  def compute_particle_colors(frame, data):
      colors = data.particles_.create_property('Color')
      colors[:,0] = data.particles.positions[:,0] / data.cell[0,0]
      colors[:,1] = 0.0
      colors[:,2] = 0.5
  pipeline.modifiers.append(compute_particle_colors)
```"""

    class NeighborMode(enum.Enum):
        """"""
        Cutoff = enum.auto()
        Bonded = enum.auto()
    neighbor_mode: ComputePropertyModifier.NeighborMode = NeighborMode.Cutoff
    'Selects the criterion for visiting neighboring particles. Supported modes are:\n\n  * ``ComputePropertyModifier.NeighborMode.Cutoff``\n  * ``ComputePropertyModifier.NeighborMode.Bonded``\n\nDefault: ``ComputePropertyModifier.NeighborMode.Cutoff``'
    cutoff_radius: float = 3.0
    "The cutoff radius up to which neighboring particles are visited to compute :py:attr:`neighbor_expressions`.\nThis parameter is only used if :py:attr:`operate_on` is set to ``'particles'``, :py:attr:`neighbor_mode` is set to ``Cutoff``,\nand :py:attr:`neighbor_expressions` has been specified.\n\nDefault: 3.0"
    expressions: Union[str, Sequence[str], dict[str, str]] = '0'
    'expressions() -> str | collections.abc.Sequence[str] | dict[str, str]\n\nThe math expression(s) to compute for each element, one for each vector component of the output property. \n\nIf the computed :py:attr:`output_property` is scalar, the expression should be a single Python string:\n\n```python\n  ComputePropertyModifier(\n      output_property = \'MyScalarProperty\',\n      expressions = \'sqrt(Position.X^2 + Position.Y^2)\'\n  )\n```\n\nIf the output property is one of the standard vector properties, which have predefined components, a sequence of expression strings should be specified, one for each component:\n\n```python\n  ComputePropertyModifier(\n      output_property = \'Color\',\n      expressions = (\'ParticleIndex/N\', \'0\', \'0.5\') # R, G, B\n  )\n```\n\nWhen creating a user-defined vector property, a dictionary of key-value pairs should be specified to set the property\'s :py:attr:`component_names`:\n\n```python\n  ComputePropertyModifier(\n      output_property = \'MyVectorProperty\',\n      expressions = {\n          \'X\': \'Position.X - Initial.X\',\n          \'Y\': \'Position.Y - Initial.Y\',\n          \'Z\': \'Position.Z - Initial.Z\'\n      }\n  )\n```\n\nFor a description of the math expression syntax, see the corresponding user manual page. \n\n.. attention::\n\n  The internal math expression parser can only handle variable identifiers made of alphanumeric characters and underscores.   OVITO property and attribute names, however, can contain arbitrary characters (except dots). Within a math expression,   properties and attributes must thus be referenced by their *mangled* names, with all spaces removed and other invalid characters replaced by   underscores. For example, if the *real* name of an input property is "*c_abc[1] Average*",   then it must be referenced in the math expression as ``c_abc_1_Average``. Furthermore, if the first character of a   name happens to be a digit, an extra underscore is prepended to the mangled name. \n\nThe modifier automatically creates a :py:class:`VectorVis` visual element for a user-defined particle property if it has 3 components that are named ``\'X\'``, ``\'Y\'``, ``\'Z\'``. This :py:class:`VectorVis` element, which is attached to the generated :py:class:`Property` object, is disabled by default and must be enabled and configured as follows if an arrow glyph visualization is desired:\n\n```python\n  pipeline.modifiers.append(ComputePropertyModifier(\n      output_property = \'InPlaneDir\',\n      expressions = {\'X\': \'cos(Angle)\', \'Y\': \'sin(Angle)\', \'Z\': \'0\'}\n  ))\n  vis = pipeline.compute().particles[\'InPlaneDir\'].vis\n  vis.enabled = True\n  vis.color = (1,0,0)\n  vis.width = 0.2\n```\n\nDefault: ``\'0\'``\n\n   Support for user-defined vector properties and automatic creation of a :py:class:`VectorVis` visual element.'
    neighbor_expressions: Union[str, Sequence[str]] = tuple()
    "The math expression(s) for the per-neighbor term(s), one for each vector component of the output particle property.\nThe number of expressions in the list must match the number of vector components of the output property.\nIf the output property is a scalar, only one expression string is required.\n\nThe neighbor expressions are evaluated for every neighboring particle and the obtained values are\nadded to the output value of the central particle. Which neighboring particles are considered\ndepends on the :py:attr:`neighbor_mode` and :py:attr:`cutoff_radius`. See :ref:`manual:particles.modifiers.compute_property.neighbor_expr`\nfor more information.\n\nNeighbor expressions are only used if :py:attr:`operate_on` is set to ``'particles'``.\n\nExample: Compute mean velocity vector, averaged over all neighbors of a particle within a 3.8 Angstrom radius:\n\n```python\n  ComputePropertyModifier(\n      output_property = 'Averaged Velocity',\n      expressions = {\n          'X': 'Velocity.X / (NumNeighbors+1)',\n          'Y': 'Velocity.Y / (NumNeighbors+1)',\n          'Z': 'Velocity.Z / (NumNeighbors+1)'\n      },\n      neighbor_expressions = (\n          'Velocity.X / (NumNeighbors+1)',\n          'Velocity.Y / (NumNeighbors+1)',\n          'Velocity.Z / (NumNeighbors+1)'\n      ),\n      cutoff_radius = 3.8\n  )\n```\n\nDefault: ``()``"
    only_selected: bool = False
    'only_selected() -> bool\n\nIf ``True``, the property is only computed for currently selected elements. In this case, the property values of unselected elements will be preserved if the output property already exists. \n\nDefault: ``False``'
    operate_on: str = 'particles'
    "operate_on() -> str\n\nSelects the kind of data elements this modifier will operate on. Valid values are: \n\n  * ``'particles'`` - to add a property to the :py:class:`Particles` container\n  * ``'bonds'`` - to add a property to the :py:class:`Bonds` container\n  * ``'voxels:<ID>'`` - to add a property to the :py:class:`VoxelGrid` container with the given :py:attr:`identifier`\n  * ``'lines:<ID>'`` - to add a property to a :py:class:`Lines` container with the given :py:attr:`identifier`\n  * ``'vectors:<ID>'`` - to add a property to a :py:class:`Vectors` container with the given :py:attr:`identifier`\n\n\nDefault: ``'particles'``"
    output_property: str = 'My property'
    "output_property() -> str\n\nThe name of the output property which the computed values will be assigned to. If the property already exists in the input data, its values will be overwritten. If the name matches one of the predefined standard properties of OVITO, it will use the prescribed data type and component names. Otherwise, a new user-defined property with data type ``float64`` is created. \n\n* List of standard particle properties\n* List of standard bond properties\n* List of standard voxel properties\n* List of standard lines properties\n* List of standard vectors properties\n\n\nDefault: ``'My property'``"

@dataclass(kw_only=True)
class ConstructSurfaceModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

Constructs the geometric surface of the solid region formed by point-like particles. The modifier generates
a :py:class:`SurfaceMesh`, which is a closed manifold consisting of triangular elements. It furthermore can compute the
surface area and the volume of individual regions enclosed by the surface mesh as well as empty pores in the interior of a structure.
See the corresponding user manual page for more information.

Basic usage example:

```python
  pipeline.modifiers.append(ConstructSurfaceModifier(radius = 2.9))
  data = pipeline.compute()
  surface_mesh = data.surfaces['surface']
```

See the :py:class:`SurfaceMesh` class for further information on how to access the mesh's vertices, faces, and spatial regions.

Construction methods:

The modifier supports two different surface construction methods, which are selected via the :py:attr:`method`
parameter. The ``AlphaShape`` method uses a Delaunay tessellation constructed on the basis of the input particle centers
to identify regions of space that cannot fully accommodate a virtual probe sphere.
These Delaunay elements are classified as being part of the filled (or solid) region.

The :py:attr:`radius` parameter controls the size of the virtual probe sphere and determines how
much detail of the morphology are resolved during the surface construction. A larger radius leads to a surface with fewer details,
reflecting only coarse features of the surface morphology. A small radius, on the other hand, will resolve finer surface features and
small pores in the interior of a solid, for example.

The second surface construction method supported by the modifier, ``GaussianDensity``, is based on the smeared-out representation of the point-like
particles as overlapping Gaussian distribution functions centered on each particle site.
The surface mesh is constructed as an isosurface of the Gaussian density field, with the threshold :py:attr:`isolevel` chosen such that
the resulting manifold roughly represents the input particle spheres. This code demonstrates how to use this second construction method:

```python
  pipeline.modifiers.append(ConstructSurfaceModifier(
      method = ConstructSurfaceModifier.Method.GaussianDensity,
      radius_scaling = 1.2,
      isolevel = 0.5))
```

Keep in mind that the result will depend on the current display radii of the particles when using the ``GaussianDensity`` method.
The effective radii can generally depend on the standard radius set in the :py:class:`ParticlesVis` element, the :py:attr:`radius` set for each
:py:class:`ParticleType`, or the per-particle property ``Radius`` if it is defined.

Volumetric region analysis:

Both methods provide the option to compute volumes and surface areas of the volumetric regions.
This allows for the identification -- and quantitative characterization -- of individual cavities or clusters formed by particles.
To enable this extra calculation step, set :py:attr:`identify_regions` to ``True``. The modifier will also output various
aggregate values for the entire system such as the total surface area and the porosity as global attributes,
see the table below.

.. attention::

  For simulation cells with one or more non-periodic boundaries, the volume of exterior regions is not
  well defined. Therefore, OVITO will report the volumes of theses regions as infinity (*inf*). Global attributes
  derived from these volumes will also be reported as either infinity or not a number (*nan*).

```python
  pipeline.modifiers.append(ConstructSurfaceModifier(
      method = ConstructSurfaceModifier.Method.AlphaShape,
      radius = 2.9,
      identify_regions = True))
  
  data = pipeline.compute()
  
  print("Aggregate quantities:")
  print(f"Surface area: {data.attributes['ConstructSurfaceMesh.surface_area']}")
  print(f"Solid volume: {data.attributes['ConstructSurfaceMesh.filled_volume']}")
  
  print("Volumetric regions:")
  regions = data.surfaces['surface'].regions
  for vol, area, filled, exterior in zip(regions['Volume'], regions['Surface Area'], regions['Filled'], regions['Exterior']):
      print(f"  volume={vol}, surface area={area}, filled={filled}, exterior={exterior}")
```

The ``AlphaShape`` method furthermore provides the option to determine for each input particle which spatial region(s) it belongs to.
For more information, see the :py:attr:`map_particles_to_regions` parameter.

Outputs:

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Surface
      -
    * - :py:class:`SurfaceMesh`
      - The surface mesh computed by the modifier. You can access it through the :py:attr:`surfaces`
        dictionary in the output :py:class:`DataCollection` under the lookup key ``"surface"``, see the
        code example below.

The following output attributes are computed only if :py:attr:`identify_regions` is turned on:

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Global attributes
      -
    * - ``ConstructSurfaceMesh.surface_area``
      - The total surface area in squared simulation units of length.
        Includes only interfaces between empty and filled regions if they modifier identifies internal boundaries (i.e. filled-filled interfaces).
    * - ``ConstructSurfaceMesh.specific_surface_area``
      - The total surface area divided by the simulation cell volume (in reciprocal units of length).
    * - ``ConstructSurfaceMesh.filled_volume``
      - The total volume of the filled region(s) bounded by the surface in cubic simulation units of length.
    * - ``ConstructSurfaceMesh.empty_volume``
      - The total volume of the empty regions in cubic simulation units of length. Includes interior pores as well as the overlap of the infinite exterior region
        with the simulation box volume in case of a finite system with open boundary conditions.
    * - ``ConstructSurfaceMesh.void_volume``
      - The total volume of the internal empty regions in cubic simulation units of length. Only includes interior pores that are *not* adjacent to an open simulation boundary.
    * - ``ConstructSurfaceMesh.filled_fraction``
      - Total volume of filled regions divided by the simulation box volume (unitless).
    * - ``ConstructSurfaceMesh.empty_fraction``
      - Total volume of empty regions divided by the simulation box volume (unitless).
    * - ``ConstructSurfaceMesh.filled_region_count``
      - Integer number of disconnected volumetric regions filled with particles.
    * - ``ConstructSurfaceMesh.empty_region_count``
      - Integer number of disconnected empty regions (any).
    * - ``ConstructSurfaceMesh.void_region_count``
      - Integer number of internal disconnected empty regions, which are not adjacent to an open simulation box boundary.
    * - ``ConstructSurfaceMesh.cell_volume``
      - The total volume of the simulation box in cubic simulation units of length. Equal to `SimulationCell.volume`.
    * - ``ConstructSurfaceMesh.region_memberships``
      - A list of arrays, specifying for each spatial region which input particles belong to it.
        This attribute is only generated if the :py:attr:`map_particles_to_regions` option has been enabled.

The following particle properties may be computed by the modifier:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Selection``
      - Will be set to 1 by the modifier for particles right at the surface of a filled region, and 0 for interior particles
        or isolated particles not forming any solid. Computed only if the :py:attr:`select_surface_particles` option is turned on.
    * - ``Surface Distance``
      - The computed distance of each particle from the closest point on the geometric surface.
        This property is only computed if the :py:attr:`compute_distances` option is turned on.
    * - ``Region``
      - Index of the spatial region (see `SurfaceMesh.regions`) each particle is located in.
        This property is only computed if the :py:attr:`map_particles_to_regions` option is turned on."""

    class Method(enum.Enum):
        """"""
        AlphaShape = enum.auto()
        GaussianDensity = enum.auto()
    compute_distances: bool = False
    "compute_distances() -> bool\n\nThis option activates the calculation of distances of the particles from the constructed surface. The computed distance of each particle is measured to the closest point on the surface mesh. The modifier will output the computed distance values as a new particle property named ``Surface Distance``. \n\nNote that the computation of distances for all particles is a very expensive operation, which can take a long time for systems with many particles or a complex surface. To select surface particles, i.e. those particles which are located exactly on the surface, it may be more efficient to use the modifier's :py:attr:`select_surface_particles` option instead. \n\nDefault: ``False``"
    grid_resolution: int = 50
    'grid_resolution() -> int\n\nSpecifies the number of grid cells along the longest dimension of the simulation cell when generating the Gaussian density grid. This parameter thus also controls the level of detail of the of the final surface mesh, which is constructed as an isosurface from the density grid data. \n\nThe parameter is only used if ``GaussianDensity`` is selected as construction :py:attr:`method`.\n\nDefault: ``50``'
    identify_regions: bool = False
    'identify_regions() -> bool\n\nThis option lets the modifier identify individual volumetric regions (filled with particles or empty) and compute their volumes and surface areas. \n\nDefault: ``False``'
    isolevel: float = 0.6
    'isolevel() -> float\n\nThe threshold value used for constructing the isosurface from the Gaussian density field. This parameter is only used if the selected construction :py:attr:`method` is set to ``GaussianDensity``. \n\nDefault: ``0.6``'
    method: ConstructSurfaceModifier.Method = Method.AlphaShape
    'method() -> ConstructSurfaceModifier.Method\n\nSelects the algorithm for constructing the surface from the input particles. The following methods are available:\n\n\n * ``ConstructSurfaceModifier.Method.AlphaShape``\n * ``ConstructSurfaceModifier.Method.GaussianDensity``\n\nSee also :ref:`manual:particles.modifiers.construct_surface_mesh.method_comparison`.\n\nIf the alpha-shape method is selected, you should also set the :py:attr:`radius` parameter to specify the size of the virtual probe sphere, which determines the level of detail of the resulting surface. \n\nDefault: ``ConstructSurfaceModifier.Method.AlphaShape``'
    only_selected: bool = False
    'only_selected() -> bool\n\nIf ``True``, the modifier acts only on selected particles and ignores other particles; if ``False``, the modifier constructs the surface around all particles.\n\nDefault: ``False``'
    radius: float = 4.0
    'radius() -> float\n\nThe radius of the virtual probe sphere used in the alpha-shape surface construction algorithm. This parameter is only used by the modifier if :py:attr:`method` is set to the default mode ``AlphaShape``. \n\nA rule of thumb is that the radius parameter should be slightly larger than the typical distance between nearest neighbor particles. \n\nDefault: ``4.0``'
    radius_scaling: float = 1.0
    'radius_scaling() -> float\n\nScaling factor applied to the input particle radii when constructing the Gaussian density field for surface generation. This parameter is only used if the selected construction :py:attr:`method` is set to ``GaussianDensity``. The optional scaling serves as a way to increase the spatial extent of the Gaussian function centered on each particle site and to increase the overlap between the Gaussians, yielding a more connected isosurface. \n\nDefault: ``1.0``'
    select_surface_particles: bool = False
    'select_surface_particles() -> bool\n\nThis option lets the modifier select those particles that are part of the constructed geometric surface. This provides a simply way of identifying surface atoms or particles. If the flag is set to ``True``, the modifier will create the ``Selection`` particle property, assigning a value of 1 to surface particles and 0 to bulk particles. This particle selection may then be used in subsequent operations in the data pipeline. \n\nSelecting surface particles is only supported by the ``AlphaShape`` construction :py:attr:`method`. For other construction methods, the setting is ignored and the modifier does not create a particle selection. \n\nAn alternative way of selecting particles that are located near the surface is to use the :py:attr:`compute_distances` option of the modifier. While the :py:attr:`select_surface_particles` option lets you identify particles located exactly on the geometric surface, the :py:attr:`compute_distances` option lets you define a distance threshold, selecting also particles slightly away from the surface. \n\nDefault: ``False``'
    smoothing_level: int = 8
    'smoothing_level() -> int\n\nThe number of times the smoothing procedure is applied to the generated surface mesh. This parameter is only used by the modifier if :py:attr:`method` is set to the default mode ``AlphaShape``. \n\nNote that the smoothing level does only affect the computed surface area but not the solid volume. That is because the solid volume is computed before smoothing the mesh. (Smoothing is supposed to be volume preserving.)\n\nDefault: ``8``'
    transfer_properties: bool = False
    'transfer_properties() -> bool\n\nThis option lets the modifier copy the property values of the particles over to the vertices of the generated surface mesh. \n\nNote: If the Gaussian density method is used, only particle properties of data type ``Float`` will be transferred to the surface. Integer properties will be skipped, because the algorithm has to blend the property values of several particles in order to compute the value at each surface vertex. In case of the alpha-shape method, all properties are transferred, including integer properties, because there is a one-to-one mapping between input particles and output mesh vertices. \n\nDefault: ``False``'
    map_particles_to_regions: bool = False
    'map_particles_to_regions() -> bool\n\nThis option lets the modifier determine for each input particles which spatial region(s) of\nthe output :py:class:`SurfaceMesh` it is part of.\nThe ``AlphaShape`` surface construction :py:attr:`method` must be selected and the :py:attr:`identify_regions` option turned on for this\nto work; otherwise the :py:attr:`map_particles_to_regions` option is ignored.\n\nTwo outputs are produced by this option:\n\n1. A per-particle property named ``Region`` is generated, which specifies for each particle which spatial region it is located in (0-based region index).\n   If a particle is located exactly on the boundary between a filled and an empty region, it will always be attributed to the filled region.\n   However, a particle located fully within an empty region will be attributed to the empty region (which can happen for solitary particles that do\n   not constitute a filled region).\n\n   ```python\n  # Construct surface mesh around the subset of currently selected particles.\n  pipeline.modifiers.append(ConstructSurfaceModifier(radius = 5.0, smoothing_level = 0,\n      only_selected = True, identify_regions = True, map_particles_to_regions = True))\n  data = pipeline.compute()\n  \n  # Count particles located inside solid regions of the surface mesh (region property "Filled" is 1).\n  particle_regions = data.particles[\'Region\']\n  regions_table = data.surfaces[\'surface\'].regions\n  n = numpy.count_nonzero(regions_table[\'Filled\'][particle_regions])\n```\n\n2. In general, a particle can belong to several regions simultaneously if it is located right on the boundary between two adjacent regions.\n   It may even be affiliated with more than one *filled* region if the surface mesh contains internal interfaces because it was constructed\n   from a subdivided input structure (see :py:class:`GrainSegmentationModifier`).\n\n   In such cases, the per-particle property ``Region`` described above is insufficient to represent the multiple region affiliations a particle can have.\n   Therefore, the modifier also generates a membership list for each spatial region, which lists the 0-based indices of all particles affiliated with that region.\n   This list of lists is stored in the output data collection as the global attribute ``ConstructSurfaceMesh.region_memberships``:\n\n   ```python\n  # Print the list of particles affiliated with each region:\n  memberships = data.attributes[\'ConstructSurfaceMesh.region_memberships\']\n  for idx, plist in enumerate(memberships):\n      print(f"Region {idx} contains {len(plist)} particles: {plist[...]}")\n```\n\nNote that empty regions will typically *not* have an empty particle membership list, because the surface particles at the interface between a filled and an empty spatial\n   region get affiliated with *both* of these regions.\n\n   You can select all particles that belong to a particular spatial region as follows:\n\n   ```python\n  region_index_to_select = 0\n  selection = data.particles_.create_property(\'Selection\', data=0)\n  selection[memberships[region_index_to_select]] = 1\n```\n\nThe mapping of particles to regions happens *before* the surface smoothing stage of the alpha-shape construction algorithm,\n    which means particles may appear slightly outside of the region they were attributed to, because the vertices of the final\n    surface mesh are displaced by the smoothing procedure. You can turn off surface smoothing by setting :py:attr:`smoothing_level` to 0.\n\n    The modifier assigns all particles to some spatial region(s) - even the unselected particles if the :py:attr:`only_selected` option is turned on\n    and the surface is being constructed only from a subset of the input particles.\n\nDefault: ``False``'
    vis: ovito.vis.SurfaceMeshVis = ovito.vis.SurfaceMeshVis()
    'vis() -> ovito.vis.SurfaceMeshVis\n\nThe :py:class:`SurfaceMeshVis` element controlling the visual appearance of the generated surface in rendered images and animations.'

@dataclass(kw_only=True)
class CoordinationAnalysisModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

Computes coordination number of each particle and the radial distribution function (RDF) for the entire system.
See the corresponding user manual page
for more information.

Inputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Position``
      - The input coordinates of the particles.
    * - ``Particle Type``
      - Required if :py:attr:`partial` is set to ``True``.
    * - ``Selection``
      - The selection state of the input particles. Only needed if :py:attr:`only_selected` is set to ``True``.

Outputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Coordination``
      - The number of neighbors of each particle.

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - :py:class:`DataTable`
      -
    * - ``coordination-rdf``
      - The RDF computed by the modifier. You can retrieve the RDF data through the :py:attr:`tables`
        dictionary of the output :py:class:`DataCollection` under the lookup key ``"coordination-rdf"``, see the
        code examples below.

Examples:

The following batch script demonstrates how to load a particle configuration, compute the RDF using the modifier
and export the data to a text file:

```python
  from ovito.io import import_file
  from ovito.modifiers import CoordinationAnalysisModifier
  
  # Load a particle dataset, apply the modifier, and evaluate pipeline.
  pipeline = import_file("input/simulation.dump")
  modifier = CoordinationAnalysisModifier(cutoff = 5.0, number_of_bins = 200)
  pipeline.modifiers.append(modifier)
  data = pipeline.compute()
  
  # Print the computed g(r) function values.
  print(data.tables['coordination-rdf'].xy())
```

The following code demonstrates how to use the :py:class:`TimeAveragingModifier` to compute the RDF for every frame of an MD simulation and 
generate a time-averaged RDF histogram. Finally, the RDF histogram is written to an output file.

```python
  from ovito.io import import_file, export_file
  from ovito.modifiers import CoordinationAnalysisModifier, TimeAveragingModifier
  import numpy
  
  # Load a simulation trajectory consisting of several frames:
  pipeline = import_file("input/simulation.dump")
  print("Number of MD frames:", pipeline.num_frames)
  
  # Insert the RDF calculation modifier into the pipeline:
  pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff = 5.0, number_of_bins = 200))
  
  # Insert the time-averaging modifier into the pipeline, which accumulates
  # the instantaneous DataTable produced by the previous modifier and computes a mean histogram.
  pipeline.modifiers.append(TimeAveragingModifier(operate_on='table:coordination-rdf'))
  
  # Data export method 1: Convert to NumPy array and write data to a text file:
  total_rdf = pipeline.compute().tables['coordination-rdf[average]'].xy()
  numpy.savetxt("output/rdf.txt", total_rdf)
  
  # Data export method 2: Use OVITO's own export function for DataTable objects:
  export_file(pipeline, "output/rdf.txt", "txt/table", key="coordination-rdf[average]")
```"""
    cutoff: float = 3.2
    'cutoff() -> float\n\nSpecifies the cutoff distance for the coordination number calculation and also the range up to which the modifier calculates the RDF. \n\nDefault: ``3.2``'
    number_of_bins: int = 200
    'number_of_bins() -> int\n\nThe number of histogram bins to use when computing the RDF.\n\nDefault: ``200``'
    only_selected: bool = False
    'only_selected() -> bool\n\nRestricts the calculation to currently selected particles. Unselected particles will be treated as if they did not exist and their ``Coordination`` value is set to zero when this option is enabled. \n\nDefault: ``False``'
    partial: bool = False
    'partial() -> bool\n\nThis modifier option requests the calculation of element-specific (partial) RDFs. The resulting RDF table will contain one tabulated g(r) function for each pair-wise combination of particle types in the input. This code example demonstrates how to access the partial RDFs computed by the modifier: \n\n```python\n  from ovito.io import import_file\n  from ovito.modifiers import CoordinationAnalysisModifier\n  \n  # Load input data.\n  pipeline = import_file("input/simulation.dump")\n  \n  # Print the list of input particle types.\n  # They are represented by ParticleType objects attached to the \'Particle Type\' particle property.\n  for t in pipeline.compute().particles.particle_types.types:\n      print("Type %i: %s" % (t.id, t.name))\n  \n  # Calculate partial RDFs:\n  pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff=5.0, number_of_bins=100, partial=True))\n  \n  # Access the output DataTable:\n  rdf_table = pipeline.compute().tables[\'coordination-rdf\']\n  \n  # The y-property of the data points of the DataTable is now a vectorial property.\n  # Each vector component represents one partial RDF.\n  rdf_names = rdf_table.y.component_names\n  \n  # Print a list of partial g(r) functions.\n  for component, name in enumerate(rdf_names):\n      print("g(r) for pair-wise type combination %s:" % name)\n      print(rdf_table.y[:,component])\n  \n  # The DataTable.xy() method yields everthing as one combined NumPy table.\n  # This includes the \'r\' values in the first array column, followed by the\n  # tabulated g(r) partial functions in the remaining columns. \n  print(rdf_table.xy())\n```\n\nDefault: ``False``'

@dataclass(kw_only=True)
class CoordinationPolyhedraModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

Constructs coordination polyhedra around selected particles for visualization purposes.
See the corresponding user manual page
for more information. A coordination polyhedron is the convex hull constructed from the neighbor atoms
of some central atom.

In order to tell the modifier which input particles should be surrounded by a coordination polyhedron,
the central particles must be selected first. The particle selection can be created by inserting a
:py:class:`SelectTypeModifier` into the data pipeline prior to this modifier.

The modifier furthermore requires input bonds connecting a central particle with its neighbors
to define which atoms the algorithm should use in the construction of the convex hulls.
You can insert the :py:class:`CreateBondsModifier` prior to this modifier to let OVITO dynamically generate
neighbor bonds between particles.

Example:

```python
  # Import simulation dataset and add it to the scene:
  pipeline = import_file("input/simulation.0.dump")
  pipeline.add_to_scene()
  
  # Select all atoms of species type 1. They will form the centers of the polyhedra.
  pipeline.modifiers.append(SelectTypeModifier(types={1}))
  # Create bonds between nearby atoms.
  pipeline.modifiers.append(CreateBondsModifier(cutoff=3.0))
  # Let the modifier construct the coordination polyhedra around selected atoms.
  modifier = CoordinationPolyhedraModifier()
  pipeline.modifiers.append(modifier)
  
  # Optional: Configure visual appearance of the polyhedra and make them semi-transparent.
  modifier.vis.surface_transparency = 0.4
  modifier.vis.highlight_edges = True
```

Modifier inputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Selection``
      - Determines the subset of particles for which coordination polyhedra should be constructed.
        You can select all particles of certain chemical type(s) by first inserting a :py:class:`SelectTypeModifier` into the pipeline.

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Bond properties
      -
    * - ``Topology``
      - The bonds list, which is used to determine the set of bonded neighbor particles serving as vertices
        for the construction of a coordination polyhedron around some central particle.
        You can let OVITO create the bond topology by inserting a :py:class:`CreateBondsModifier` into the pipeline first.

Modifier outputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Polyhedral mesh
      -
    * - ``coord-polyhedra``
      - The polyhedral :py:class:`SurfaceMesh` generated by the modifier. You can retrieve it from the :py:attr:`surfaces`
        dictionary of the :py:class:`DataCollection`."""
    vis: ovito.vis.SurfaceMeshVis = ovito.vis.SurfaceMeshVis(show_cap=False, smooth_shading=False, surface_transparency=0.25, title='Polyhedra')
    'vis() -> ovito.vis.SurfaceMeshVis\n\nThe :py:class:`SurfaceMeshVis` element controlling the visual appearance of the generated polyhedra in rendered images and animations.'

@dataclass(kw_only=True)
class CreateBondsModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

Creates bonds between nearby particles.
See the corresponding user manual page
for more information.

Inputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Position``
      - The xyz coordinates of the input particles.
    * - ``Particle Type``
      - The particle type information, which is used only if :py:attr:`mode` is set to ``VdWRadius`` or ``Pairwise``.
    * - ``Molecule Identifier``
      - The assignment of atoms to molecules, which is used only if :py:attr:`intra_molecule_only` is set to ``True``.

Outputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Bond properties
      -
    * - ``Topology``
      - The modifier will create new bond topology entries and append them to the property arrays in an existing :py:class:`Bonds`
        object; or it creates a new :py:class:`Bonds` instance if necessary.
    * - ``Periodic Image``
      - Stores the transitions of each bond through the faces of a periodic simulation cell if the bond connects
        two particles from different periodic images of the system.
    * - ``Bond Type``
      - The type ID information that is assigned to newly created bonds according to the modifier's :py:attr:`bond_type` field.

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Global attributes
      -
    * - ``CreateBonds.num_bonds``
      - The number of bonds that exists after the modifier's operation."""

    class Mode(enum.Enum):
        """"""
        Uniform = enum.auto()
        VdWRadius = enum.auto()
        Pairwise = enum.auto()
    bond_type: ovito.data.BondType = ovito.data.BondType()
    'bond_type() -> ovito.data.BondType\n\nThe :py:class:`BondType` that will be assigned to the newly created bonds. This lets you control the display color of the new bonds.'
    cutoff: float = 3.2
    'cutoff() -> float\n\nThe upper cutoff distance for the creation of bonds between particles. This parameter is only used if :py:attr:`mode` is ``Uniform``. \n\nDefault: ``3.2``'
    intra_molecule_only: bool = False
    'intra_molecule_only() -> bool\n\nIf this option is set to true, the modifier will create bonds only between atoms that belong to the same molecule (i.e. which have the same molecule ID assigned to them).\n\nDefault: ``False``'
    lower_cutoff: float = 0.0
    'lower_cutoff() -> float\n\nThe minimum bond length. No bonds will be created between particles whose distance is below this threshold. \n\nDefault: ``0.0``'
    mode: CreateBondsModifier.Mode = Mode.Uniform
    'mode() -> CreateBondsModifier.Mode\n\nControls the cutoff criterion for creating bonds. Valid cutoff modes are:\n\n  * ``CreateBondsModifier.Mode.Uniform``\n  * ``CreateBondsModifier.Mode.VdWRadius``\n  * ``CreateBondsModifier.Mode.Pairwise``\n\n\nMode ``Uniform`` uses a single uniform :py:attr:`cutoff` distance for creating bonds, which is independent of the types of the particles. Mode ``VdWRadius`` uses a distance cutoff that is derived from the :py:attr:`vdw_radius` (Van der Waals radius) of the :py:class:`ParticleType` of the two particles involved. Mode ``Pairwise`` lets you specify a separate cutoff distance for each pairwise combination of particle types using the :py:meth:`set_pairwise_cutoff` method. \n\nDefault: ``CreateBondsModifier.Mode.Uniform``'
    prevent_hh_bonds: bool = True
    "prevent_hh_bonds() -> bool\n\nControls whether the modifier should *not* generate bonds between pairs of hydrogen atoms. This flag only applies if :py:attr:`.mode` is set to ``VdWRadius`` and the van der Waals radii of the atom types are used for generating pair-wise bonds. A particle is considered a hydrogen atom if its :py:class:`ParticleType`'s name is ``'H'``. \n\nDefault: ``True``"
    vis: ovito.vis.BondsVis = ovito.vis.BondsVis()
    'vis() -> ovito.vis.BondsVis\n\nThe :py:class:`BondsVis` object controlling the visual appearance of the bonds created by this modifier.'

    def get_pairwise_cutoff(self, type_a: Union[int, str], type_b: Union[int, str]) -> float:
        """Returns the pair-wise cutoff distance that was previously set for a specific pair of particle types using the :py:meth:`set_pairwise_cutoff` method. 

:param type_a: The :py:attr:`name` or numeric :py:attr:`id` of the first particle type
:param type_b: The :py:attr:`name` or numeric :py:attr:`id` of the second particle type
:return: The cutoff distance set for the type pair. Returns zero if no cutoff has been set for the pair."""
        ...

    def set_pairwise_cutoff(self, type_a: Union[int, str], type_b: Union[int, str], cutoff: float) -> None:
        """Sets the cutoff range for creating bonds between a specific pair of particle types. This information is only used if :py:attr:`mode` is set to ``Pairwise``.

:param type_a: The :py:attr:`name` or numeric :py:attr:`id` of the first particle type
:param type_b: The :py:attr:`name` or numeric :py:attr:`id` of the second particle type
:param cutoff: The cutoff distance to be used by the modifier for the type pair


If you want no bonds to be created between a pair of types, set the corresponding cutoff radius to zero (which is the default)."""
        ...

@dataclass(kw_only=True)
class CreateIsosurfaceModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

Constructs an isosurface for a scalar field defined on a three-dimensional :py:class:`VoxelGrid`.
See the corresponding user manual page
for more information on this modifier.

The modifier takes an existing :py:class:`VoxelGrid` as input, for example a charge density
field loaded from a simulation file or a three-dimensional voxel grid dynamically
computed by the :py:class:`SpatialBinningModifier`. The input voxel grid to be used in isosurface computation
must be selected by setting the :py:attr:`operate_on` field to the string ``'voxels:'`` followed by the
:py:attr:`identifier` of the input voxel grid. Furthermore, you need to specify
the name of the field :py:class:`Property` for which to construct the isosurface, because
multiple fields may be defined on the same voxel grid.

    ```python
  # Import a charge density field from a VASP calculation:
  pipeline = import_file("input/CHGCAR.nospin.gz")
  pipeline.add_to_scene()
  
  # The identifier of the imported VoxelGrid is 'charge-density',
  # and the field property defined on the grid is named 'Charge Density'.
  print(pipeline.compute().grids['charge-density'])
  
  # Set up the isosurface modifier:
  modifier = CreateIsosurfaceModifier(
      operate_on = 'voxels:charge-density',
      property = 'Charge Density',
      isolevel = 0.00014)
  pipeline.modifiers.append(modifier)
  
  # Adjust visual appearance of the isosurface in rendered images:
  modifier.vis.show_cap = False
  modifier.vis.surface_transparency = 0.4
```

The following example demonstrates how to construct an isosurface for a field dynamically computed by
the :py:class:`SpatialBinningModifier`:

    ```python
  pipeline.modifiers.append(SpatialBinningModifier(
      property = 'Particle Type',
      direction = SpatialBinningModifier.Direction.XYZ,
      bin_count = (30, 30, 30),
      reduction_operation = SpatialBinningModifier.Operation.Mean
  ))
  pipeline.modifiers.append(CreateIsosurfaceModifier(
      operate_on = 'voxels:binning',
      property = 'Particle Type',
      isolevel = 1.5
  ))
```

Modifier outputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - :py:class:`SurfaceMesh`
      -
    * - ``isosurface``
      - The :py:class:`SurfaceMesh` constructed by the modifier. You can retrieve it from the :py:attr:`surfaces`
        dictionary of the pipeline's output :py:class:`DataCollection`.

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - :py:class:`DataTable`
      -
    * - ``isosurface-histogram``
      - A histogram of the input field values. You can retrieve this :py:class:`DataTable` from the :py:attr:`tables`
        dictionary of the pipeline's output :py:class:`DataCollection`."""
    isolevel: float = 0.0
    'isolevel() -> float\n\nThe value at which to create the isosurface.\n\nDefault: ``0.0``'
    operate_on: str = 'voxels:'
    "operate_on() -> str\n\nSpecifies the :py:class:`VoxelGrid` this modifier should operate on. Set this to the prefix string ``'voxels:'`` followed by the :py:attr:`identifier` of the grid. \n\nNote: You can use the Python statement ``print(pipeline.compute().grids)`` to find out what the identifiers of the voxel :py:attr:`grids` in your data pipeline are. \n\nDefault: ``'voxels:'``"
    property: str = ''
    "property() -> str\n\nThe name of the :py:class:`Property` in the selected input :py:class:`VoxelGrid` for which to construct the isosurface. Note that this must be a scalar property. If the input grid property is a vector property, you need to explicitly specify which vector component to use, e.g. ``'Dipole Orientation.Z'``."
    transfer_values: bool = False
    'transfer_values() -> bool\n\nThis option lets the modifier copy auxiliary properties of the voxel field (aside from the primary field quantity for which the isosurface is being constructed) over to the vertices of the generated isosurface mesh. The local field value at each isosurface vertex will be computed from the voxel values using trilinear interpolation. Subsequently, the pseudo-coloring capability of the :py:attr:`vis` element can be used to color the generated isosurface based on a secondary field quantity. \n\nDefault: ``False``'
    identify_regions: bool = False
    'identify_regions() -> bool\n\nThis option lets the modifier identify volumetric regions, i.e. regions of space where the field value is either below or above the iso-level, and compute their individual volumes and surface areas. These properties are calculated individually for each disconnected spatial region bounded by the isosurface. \n\nDefault: ``False``'
    smoothing_level: int = 0
    'smoothing_level() -> int\n\nThe number of times the mesh smoothing procedure is applied to the generated isosurface to even out surface steps resulting from the discrete nature of the input voxel grid. \n\nDefault: ``0``'
    vis: ovito.vis.SurfaceMeshVis = ovito.vis.SurfaceMeshVis(show_cap=False, title='Isosurface')
    'vis() -> ovito.vis.SurfaceMeshVis\n\nThe :py:class:`SurfaceMeshVis` element controlling the visual appearance of the generated isosurface in rendered images.'

@dataclass(kw_only=True)
class DeleteSelectedModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

This modifier deletes the currently selected elements.
See also the corresponding user manual page for more information.

Inputs:

The modifier can operate on any combination of the following data elements:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Data element
      -
    * - ``particles``
      - Deletes all particles with a non-zero value of the ``Selection`` property.
    * - ``bonds``
      - Deletes all bonds with a non-zero value of the ``Selection`` property.
    * - ``surface_regions``
      - Deletes all selected spatial regions (including the corresponding mesh faces) from a :py:class:`SurfaceMesh`.
    * - ``lines``
      - Deletes all vertices of a :py:class:`Lines` object having a non-zero ``Selection`` property value.
    * - ``vectors``
      - Deletes all vector glyphs of a :py:class:`Vectors` object having a non-zero ``Selection`` property value.

By default the modifier will act on all data element types simultaneously. You can restrict it to a subset by setting the :py:attr:`operate_on` field."""
    operate_on: MutableSet[str] = {'particles', 'bonds', 'surface_regions'}
    "operate_on() -> collections.abc.MutableSet[str]\n\nA set of strings specifying the kinds of data elements this modifier should operate on. By default the set contains all data element types supported by the modifier. \n\nDefault: ``{'particles', 'bonds', 'surface_regions', 'lines', 'vectors'}``"

@dataclass(kw_only=True)
class DislocationAnalysisModifier(StructureIdentificationModifier):
    """Base: :py:class:`ovito.modifiers.StructureIdentificationModifier`

This analysis modifier extracts all dislocations in a crystal and converts them to continuous line segments.
The computational method behind this is called *Dislocation Extraction Algorithm* (DXA) and is described
in the paper `[MSMSE 20 (2012), 085007] <http://stacks.iop.org/0965-0393/20/085007>`__.
See also the corresponding user manual page for this modifier.

The extracted dislocation lines are output as a :py:class:`DislocationNetwork` object by the modifier
and can be accessed through the `DataCollection.dislocations` field
after the modification pipeline has been evaluated. This is demonstrated in the example script below.

Furthermore, you can use the :py:func:`export_file` function to write the dislocation lines
to a so-called CA file. The CA file format is described on this page
of the OVITO user manual.

Inputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Position``
      - The coordinates of the input particles.
    * - ``Selection``
      - The selection state of the input particles. Only needed if :py:attr:`~StructureIdentificationModifier.only_selected` is ``True``.

Outputs:

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Dislocations
      -
    * - :py:class:`DislocationNetwork`
      - The dislocation lines found by the modifier. You can access this data object through the :py:attr:`dislocations`
        field of the output :py:class:`DataCollection`, see the
        code example below.

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Global attributes
      -
    * - ``DislocationAnalysis.total_line_length``
      - Total length of all dislocation lines found by the DXA (in simulation units).
    * - ``DislocationAnalysis.length.1/n<ijk>``
      - A set of attributes indicating the length of dislocations broken down by Burgers vector type.
        For example, the attribute ``DislocationAnalysis.length.1/6<112>`` specifies the total line length of Shockley partials found by the DXA.
    * - ``DislocationAnalysis.length.other``
      - Total length of dislocation lines with an unusual Burgers vector that do not belong to any of the predefined standard dislocation types.
    * - ``DislocationAnalysis.cell_volume``
      - The volume of the input simulation cell. You can use it to calculate the dislocation density from the line length.
    * - ``DislocationAnalysis.counts.OTHER``
      - Number of particles not matching any of the known structural types.
    * - ``DislocationAnalysis.counts.FCC``
      - Number of particles identified as face-centered cubic structure.
    * - ``DislocationAnalysis.counts.HCP``
      - Number of particles identified as hexagonal close packed structure.
    * - ``DislocationAnalysis.counts.BCC``
      - Number of particles identified as body-centered cubic structure.
    * - ``DislocationAnalysis.counts.CubicDiamond``
      - Number of particles identified as cubic diamond structure.
    * - ``DislocationAnalysis.counts.HexagonalDiamond``
      - Number of particles identified as hexagonal diamond structure.

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Structure Type``
      - The structure type assigned to each particle by the CNA sub-algorithm, encoded as an integer value:

        ============= =========================================================
        Value         Python constant
        ============= =========================================================
        0             ``DislocationAnalysisModifier.Lattice.Other``
        1             ``DislocationAnalysisModifier.Lattice.FCC``
        2             ``DislocationAnalysisModifier.Lattice.HCP``
        3             ``DislocationAnalysisModifier.Lattice.BCC``
        4             ``DislocationAnalysisModifier.Lattice.CubicDiamond``
        5             ``DislocationAnalysisModifier.Lattice.HexagonalDiamond``
        ============= =========================================================
    * - ``Color``
      - A per-particle color representing the identified structure type (only if :py:attr:`~StructureIdentificationModifier.color_by_type` is ``True``).
    * - ``Cluster``
      - The ID of the cluster each atom has been assigned to. A "cluster" is a contiguous group of atoms, all being of the same
        crystalline structure type. Non-crystalline atoms are assigned to cluster ID 0.
    * - ``Dislocation``
      - Only if :py:attr:`mark_dislocation_core_atoms` is enabled. This particle property indicates the zero-based ID of the dislocation line each atom has been associated with.
        Atoms that are not part of a dislocation core are assigned the value -1.

Example:

```python
  from ovito.io import import_file, export_file
  from ovito.modifiers import DislocationAnalysisModifier
  from ovito.data import DislocationNetwork
  
  import ovito
  ovito.enable_logging()
  
  pipeline = import_file("input/simulation.dump")
  
  # Extract dislocation lines from a crystal with diamond structure:
  modifier = DislocationAnalysisModifier()
  modifier.input_crystal_structure = DislocationAnalysisModifier.Lattice.CubicDiamond
  pipeline.modifiers.append(modifier)
  data = pipeline.compute()
  
  total_line_length = data.attributes['DislocationAnalysis.total_line_length']
  cell_volume = data.attributes['DislocationAnalysis.cell_volume']
  print("Dislocation density: %f" % (total_line_length / cell_volume))
  
  # Print list of dislocation lines:
  print("Found %i dislocation lines" % len(data.dislocations.lines))
  for line in data.dislocations.lines:
      print("Dislocation %i: length=%f, Burgers vector=%s" % (line.id, line.length, line.true_burgers_vector))
      print(line.points)
  
  # Export dislocation lines to a CA file:
  export_file(pipeline, "output/dislocations.ca", "ca")
  
  # Or export dislocations to a ParaView VTK file:
  export_file(pipeline, "output/dislocations.vtk", "vtk/disloc")
```"""

    class Lattice(enum.Enum):
        """"""
        FCC = enum.auto()
        HCP = enum.auto()
        BCC = enum.auto()
        CubicDiamond = enum.auto()
        HexagonalDiamond = enum.auto()
    circuit_stretchability: int = 9
    'circuit_stretchability() -> int\n\nThe number of steps by which a Burgers circuit can stretch while it is being advanced along a dislocation line.\n\nDefault: ``9``'
    color_by_type: bool = True
    'color_by_type() -> bool\n\nControls whether the modifier assigns a color to each particle based on the identified structure type. \n\nDefault: ``True``'
    defect_mesh_smoothing_level: int = 8
    'defect_mesh_smoothing_level() -> int\n\nSpecifies the number of iterations of the surface smoothing algorithm to perform when post-processing the extracted defect mesh.\n\nDefault: ``8``'
    defect_vis: ovito.vis.SurfaceMeshVis = ovito.vis.SurfaceMeshVis(reverse_orientation=True, cap_transparency=0.5, title='Defect mesh')
    'defect_vis() -> ovito.vis.SurfaceMeshVis\n\nThe :py:class:`SurfaceMeshVis` element controlling the visual representation of the generated defect mesh.'
    disloc_vis: ovito.vis.DislocationVis = ovito.vis.DislocationVis()
    'disloc_vis() -> ovito.vis.DislocationVis\n\nThe :py:class:`DislocationVis` element controlling the visual representation of the generated dislocation lines.'
    input_crystal_structure: DislocationAnalysisModifier.Lattice = Lattice.FCC
    'input_crystal_structure() -> DislocationAnalysisModifier.Lattice\n\nThe type of crystal to analyze. Must be one of: \n\n  * ``DislocationAnalysisModifier.Lattice.FCC``\n  * ``DislocationAnalysisModifier.Lattice.HCP``\n  * ``DislocationAnalysisModifier.Lattice.BCC``\n  * ``DislocationAnalysisModifier.Lattice.CubicDiamond``\n  * ``DislocationAnalysisModifier.Lattice.HexagonalDiamond``\n\n\nDefault: ``DislocationAnalysisModifier.Lattice.FCC``'
    line_coarsening_enabled: bool = True
    'line_coarsening_enabled() -> bool\n\nFlag that enables the coarsening of extracted dislocation lines, which reduces the number of sample points along the lines.\n\nDefault: ``True``'
    line_point_separation: float = 2.5
    'line_point_separation() -> float\n\nSets the desired distance between successive sample points along the dislocation lines, measured in multiples of the interatomic spacing. This parameter controls the amount of coarsening performed during post-processing of dislocation lines.\n\nDefault: ``2.5``'
    line_smoothing_enabled: bool = True
    'line_smoothing_enabled() -> bool\n\nFlag that enables the smoothing of extracted dislocation lines after they have been coarsened.\n\nDefault: ``True``'
    line_smoothing_level: int = 1
    'line_smoothing_level() -> int\n\nThe number of iterations of the line smoothing algorithm to perform.\n\nDefault: ``1``'
    only_perfect_dislocations: bool = False
    'only_perfect_dislocations() -> bool\n\nThis flag controls whether the algorithm should extract only perfect dislocations (and no partial dislocations, which is normally done for FCC/HCP and diamond lattices). Make sure you set the :py:attr:`circuit_stretchability` parameter to a higher value when activating this option, because larger Burgers circuits are needed to identify dissociated dislocations with wide cores. \n\nDefault: ``False``'
    mark_dislocation_core_atoms: bool = False
    "mark_dislocation_core_atoms() -> bool\n\nLets the algorithm identify and mark atoms belonging to the individual dislocation cores. A new particle property ``Dislocation`` will be created containing a zero-based dislocation identifier if the atom belongs to a dislocation's core region. Atoms not belonging to any dislocation core are assigned the special value -1.\n\nDefault: ``False``"
    only_selected: bool = False
    'only_selected() -> bool\n\nSet this to ``True`` to perform the analysis on selected particles only. Particles that are not selected will be treated as if they did not exist. Use the :py:class:`SelectTypeModifier`, for example, to restrict the crystal structure identification to a sub-lattice formed by one species of particles in a multi-component system. \n\nDefault: ``False``'
    trial_circuit_length: int = 14
    'trial_circle_length() -> int\n\nThe maximum length of trial Burgers circuits constructed by the DXA to discover dislocations. The length is specified in terms of the number of atom-to-atom steps.\n\nDefault: ``14``'

@dataclass(kw_only=True)
class ElasticStrainModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

This modifier computes the atomic-level elastic strain and deformation gradient tensors in crystalline systems.
See also the corresponding user manual page for this modifier.

The modifier first performs an identification of the local crystal structure and stores the results in the ``Structure Type`` particle
property. Possible structure type values are listed under the :py:attr:`input_crystal_structure` property.
Atoms that do not form a crystalline structure or which are part of defects are assigned the special type ``OTHER`` (=0).
For these atoms the local elastic deformation cannot be computed.

If :py:attr:`calculate_deformation_gradients` is set to true, the modifier outputs a new particle property named ``Elastic Deformation Gradient``,
which contains the per-atom elastic deformation gradient tensors. Each tensor has nine components stored in column-major order.
Atoms for which the elastic deformation gradient could not be determined (i.e. which are classified as ``OTHER``) will be assigned the null tensor.

If :py:attr:`calculate_strain_tensors` is set to true, the modifier outputs a new particle property named ``Elastic Strain``,
which contains the per-atom elastic strain tensors. Each symmetric strain tensor has six components stored in the order XX, YY, ZZ, XY, XZ, YZ.
Atoms for which the elastic strain tensor could not be determined (i.e. which are classified as ``OTHER``) will be assigned the null tensor.

Furthermore, the modifier generates a particle property ``Volumetric Strain``, which stores the trace divided by three of the local elastic strain tensor.
Atoms for which the elastic strain tensor could not be determined (i.e. which are classified as ``OTHER``) will be assigned a value of zero.

Inputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Position``
      - The coordinates of the input particles.
    * - ``Selection``
      - The selection state of the input particles. Only needed if :py:attr:`~StructureIdentificationModifier.only_selected` is ``True``.

Outputs:


.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Elastic Deformation Gradient``
      - The computed per-atom elastic deformation gradient tensors.
        Each tensor consists of 9 components stored in column-major order.
        All components will be set to zero for atoms for which no elastic deformation tensor could be determined (because they were classified as ``OTHER``).
        This output property is only generated if :py:attr:`calculate_deformation_gradients` is set to ``True``.
    * - ``Elastic Strain``
      - The computed per-atom elastic strain tensor.
        Each symmetric strain tensor has six components stored in the order XX, YY, ZZ, XY, XZ, YZ.
        All components will be set to zero for atoms for which no elastic strain tensor could be determined (because they were classified as ``OTHER``).
        This output property is only generated if :py:attr:`calculate_strain_tensors` is set to ``True``.
    * - ``Volumetric Strain``
      - Scalar particle property which is set to one third of the trace of the computed local elastic strain tensor.
        Atoms for which no elastic strain tensor could be determined (because they were classified as ``OTHER``) will
        have a volumetric strain value of zero.
    * - ``Structure Type``
      - The structure type determined by the algorithm for each particle, encoded as an integer value:

        ============= =========================================================
        Value         Python constant
        ============= =========================================================
        0             ``ElasticStrainModifier.Lattice.OTHER``
        1             ``ElasticStrainModifier.Lattice.FCC``
        2             ``ElasticStrainModifier.Lattice.HCP``
        3             ``ElasticStrainModifier.Lattice.BCC``
        4             ``ElasticStrainModifier.Lattice.CubicDiamond``
        5             ``ElasticStrainModifier.Lattice.HexagonalDiamond``
        ============= =========================================================
    * - ``Color``
      - A per-particle color representing the identified structure type (only if :py:attr:`~StructureIdentificationModifier.color_by_type` is ``True``).
    * - ``Cluster``
      - The ID of the cluster each atom has been assigned to. A "cluster" is a contiguous group of atoms, all being of the same
        crystalline structure type. Non-crystalline atoms are assigned the invalid cluster ID 0."""

    class Lattice(enum.IntEnum):
        """A simple attribute-based namespace."""
        FCC = enum.auto()
        HCP = enum.auto()
        BCC = enum.auto()
        CubicDiamond = enum.auto()
        HexagonalDiamond = enum.auto()
    axial_ratio: float = math.sqrt(8 / 3)
    'axial_ratio() -> float\n\nThe *c/a* ratio of the ideal unit cell for crystals with hexagonal symmetry.\n\nDefault: ``sqrt(8/3)``'
    calculate_deformation_gradients: bool = False
    'calculate_deformation_gradients() -> bool\n\nFlag that enables the output of the calculated elastic deformation gradient tensors. The per-particle tensors will be stored in a new particle property named ``Elastic Deformation Gradient`` with nine components (stored in column-major order). Particles for which the local elastic deformation cannot be calculated, are assigned the null tensor. \n\nDefault: ``False``'
    calculate_strain_tensors: bool = True
    'calculate_strain_tensors() -> bool\n\nFlag that enables the calculation and out of the elastic strain tensors. The symmetric strain tensors will be stored in a new particle property named ``Elastic Strain`` with six components (XX, YY, ZZ, XY, XZ, YZ). \n\nDefault: ``True``'
    input_crystal_structure: ElasticStrainModifier.Lattice = Lattice.FCC
    'input_crystal_structure() -> ElasticStrainModifier.Lattice\n\nThe type of crystal to analyze. Must be one of: \n\n  * ``ElasticStrainModifier.Lattice.FCC``\n  * ``ElasticStrainModifier.Lattice.HCP``\n  * ``ElasticStrainModifier.Lattice.BCC``\n  * ``ElasticStrainModifier.Lattice.CubicDiamond``\n  * ``ElasticStrainModifier.Lattice.HexagonalDiamond``\n\n\nDefault: ``ElasticStrainModifier.Lattice.FCC``'
    lattice_constant: float = 1.0
    'lattice_constant() -> float\n\nLattice constant (*a*:sub:`0`) of the ideal unit cell.\n\nDefault: ``1.0``'
    push_strain_tensors_forward: bool = True
    'push_strain_tensors_forward() -> bool\n\nSelects the frame in which the elastic strain tensors are calculated. \n\nIf true, the *Eulerian-Almansi* finite strain tensor is computed, which measures the elastic strain in the global coordinate system (spatial frame). \n\nIf false, the *Green-Lagrangian* strain tensor is computed, which measures the elastic strain in the local lattice coordinate system (material frame). \n\nDefault: ``True``'

@dataclass(kw_only=True)
class ExpandSelectionModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

Expands the current particle selection by selecting particles that are neighbors of already selected particles.
See the corresponding user manual page for more information.

Inputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Selection``
      - The selection state of the input particles.
    * - ``Position``
      - The coordinates of the input particles (only used if :py:attr:`mode` is ``Cutoff`` or ``Nearest``).

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Bond properties
      -
    * - ``Topology``
      - The list of bonds (only used if :py:attr:`mode` is ``Bonded``).

Outputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Selection``
      - The output particle selection.

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Global attributes
      -
    * - ``ExpandSelection.num_added``
      - The number of particles that got newly selected by the modifier."""

    class ExpansionMode(enum.Enum):
        """"""
        Cutoff = enum.auto()
        Nearest = enum.auto()
        Bonded = enum.auto()
    cutoff: float = 3.2
    'cutoff() -> float\n\nThe maximum distance up to which particles are selected around already selected particles. This parameter is only used if :py:attr:`mode` is set to ``ExpansionMode.Cutoff``.\n\nDefault: ``3.2``'
    iterations: int = 1
    'iterations() -> int\n\nControls how many iterations of the modifier are executed. This can be used to select neighbors of neighbors up to a certain recursive depth.\n\nDefault: ``1``'
    mode: ExpandSelectionModifier.ExpansionMode = ExpansionMode.Cutoff
    'mode() -> ExpandSelectionModifier.ExpansionMode\n\nSelects the mode of operation, i.e., how the modifier extends the selection around already selected particles. Valid values are:\n\n  * ``ExpandSelectionModifier.ExpansionMode.Cutoff``\n  * ``ExpandSelectionModifier.ExpansionMode.Nearest``\n  * ``ExpandSelectionModifier.ExpansionMode.Bonded``\n\n\nDefault: ``ExpandSelectionModifier.ExpansionMode.Cutoff``'
    num_neighbors: int = 1
    'num_neighbors() -> int\n\nThe number of nearest neighbors to select around each already selected particle. This parameter is only used if :py:attr:`mode` is set to ``ExpansionMode.Nearest``.\n\nDefault: ``1``'

@dataclass(kw_only=True)
class ExpressionSelectionModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

Selects data elements, e.g. particles, based on a user-defined Boolean expression.
This modifier can be inserted into a :py:class:`Pipeline` to create a selection
set. which subsequent modifiers in the pipeline can operate on.

See also the corresponding user manual page for more information on this modifier.
The modifier can select different kinds of data elements depending on the value of its :py:attr:`operate_on` property:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Operate on
      - Description
    * - ``particles`` (default)
      - Selects particles based on a user-defined criterion.
    * - ``bonds``
      - Selects bonds based on a user-defined criterion.
    * - ``surface_regions``
      - Selects volumetric regions of a :py:class:`SurfaceMesh` based on a user-defined criterion.
    * - ``lines``
      - Selects line vertices of a :py:class:`Lines` object based on a user-defined criterion.
    * - ``vectors``
      - Selects vector glyphs of a :py:class:`Vectors` object based on a user-defined criterion.

Outputs:

The modifier sets or updates the value of the ``Selection`` property of each data element: 1 if
the expression evaluates to true for an element; 0 otherwise.

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Global attributes
      -
    * - ``ExpressionSelection.count``
      - The number of data elements (particles/bonds) that have been selected by the modifier.

Usage example:

The following code example demonstrates how to select all atoms whose potential energy
exceeds a threshold value of -3.6 eV. After executing the pipeline, the number of atoms that got selected
may be queried by inspecting the global attribute ``ExpressionSelection.count``,
which was created by the modifier as summarized output:

```python
  pipeline.modifiers.append(ExpressionSelectionModifier(expression = 'PotentialEnergy > -3.6'))
  data = pipeline.compute()
  print(data.attributes['ExpressionSelection.count'])
```

Note that, in many cases, the :py:class:`ExpressionSelectionModifier` class can be replaced with an equivalent
Python modifier function, which creates the selection directly.
The example above can be reimplemented without a :py:class:`ExpressionSelectionModifier` as follows:

```python
  def select_high_energy_atoms(frame, data):
      selection = data.particles['Potential Energy'] > -3.6
      data.particles_.create_property('Selection', data=selection)
  pipeline.modifiers.append(select_high_energy_atoms)
```

In this second version, the selection criterion is specified in terms of a native Python expression operating on a NumPy array
containing particle property values, which means you are no longer limited to the simple language supported by :py:class:`ExpressionSelectionModifier`.
NumPy expressions allow you to directly and efficiently compute the ``Selection`` array and thereby create the output selection set.
See :ref:`writing_custom_modifiers` for more information.

Furthermore, if you just want to count the number of particles that fulfill a certain criterion, it may not even be necessary to create a selection.
You can directly use NumPy's `count_nonzero() <https://numpy.org/doc/stable/reference/generated/numpy.count_nonzero.html>`__
function to determine the number of particles for which some Boolean expression evaluates to true:

```python
  data = pipeline.compute()
  print(numpy.count_nonzero(data.particles['Potential Energy'] > -3.6))
```"""
    expression: str = ''
    'expression() -> str\n\nThe Boolean expression to be evaluated by the modifier for each data element, formatted as a Python string. The expression syntax is documented in :ref:`OVITO\'s user manual <manual:particles.modifiers.expression_select>`; it is *not* regular Python syntax. The string\'s contents will be parsed by OVITO\'s integrated math expression interpreter (`muparser <https://beltoforion.de/en/muparser/>`__). \n\n.. attention::\n\n  Keep in mind that the Boolean expression is subject to special syntax rules.   The internal expression parser can only handle variable names consisting of alphanumeric characters and underscores.   OVITO\'s property and attribute names, however, can contain arbitrary characters (except dots). Within a Boolean expression,   properties and attributes must thus be referenced by their *mangled* names, with all spaces removed and other invalid characters replaced by   underscores. For example, if the original name of an input property is "*c_abc[1] Average*",   then it must be referenced as ``c_abc_1_Average`` in the Boolean expression. Furthermore, if the first character of an   identifier happens to be a digit, an extra underscore is prepended to the mangled variable name. \n\nYou can incorporate the values of regular Python variables in the Boolean expression simply by using Python\'s `string formatting <https://docs.python.org/3/tutorial/inputoutput.html#formatted-string-literals>`__ technique, e.g.: \n\n```python\n  threshold = float(sys.argv[1])\n  pipeline.modifiers.append(ExpressionSelectionModifier(expression = f"ShearStrain>{threshold}"))\n```\n\nKeep in mind though that the string the modifier uses now contains a *constant literal* in place of the original Python variable. In case the Python value changes, it thus is necessary to update the expression string *after* the modifier was inserted into the pipeline. Re-evaluating the data pipeline will recompute the selection, letting you perform several calculations in a for-loop, for example: \n\n```python\n  modifier = ExpressionSelectionModifier()\n  pipeline.modifiers.append(modifier)\n  for threshold in numpy.arange(0.0, 1.0, 0.05):\n      modifier.expression = f"ShearStrain>{threshold}"\n      data = pipeline.compute()\n      print(threshold, data.attributes[\'ExpressionSelection.count\'])\n```'
    operate_on: str = 'particles'
    "operate_on() -> str\n\nSelects the kind of data elements this modifier should operate on. One of the following strings: \n\n  * ``'particles'``\n  * ``'bonds'``\n  * ``'surface_regions'``\n  * ``'lines'``\n  * ``'vectors'``\n\n\nDefault: ``'particles'``"

@dataclass(kw_only=True)
class FreezePropertyModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

This modifier copies the values of a source property from some reference animation frame (frame 0 by default) to the current animation frame.
It allows preserving a particle or bond property that would otherwise change with time.
See also the corresponding user manual page for more information.
The modifier can operate on different data elements:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Data element
      -
    * - ``particles``
      - Freezes a particle property.
    * - ``bonds``
      - Freezes a bond property.
    * - ``voxels``
      - Freezes a voxel grid property.

By default the modifier will operate on a particle property. You can change this by setting the :py:attr:`operate_on` field.

Example:

.. literalinclude:: ../example_snippets/freeze_property_modifier.py
   :emphasize-lines: 12-14"""
    destination_property: str = ''
    "destination_property() -> str\n\nThe name of the output property that should be created by the modifier. It may be the same as :py:attr:`source_property`. If the destination property already exists in the modifier's input, the values are overwritten."
    freeze_at: int = 0
    "freeze_at() -> int\n\nThe animation frame number at which to freeze the input property's values. \n\nDefault: ``0``"
    operate_on: str = 'particles'
    "operate_on() -> str\n\nSelects the kind of properties this modifier should operate on. Supported values are: ``'particles'``, ``'bonds'``, ``'voxels'``. \n\nDefault: ``'particles'``"
    source_property: str = ''
    'source_property() -> str\n\nThe name of the input property that should be evaluated by the modifier on the animation frame specified by :py:attr:`freeze_at`.'

@dataclass(kw_only=True)
class GenerateTrajectoryLinesModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

This modifier periodically samples the time-dependent positions of the particles to produce a :py:class:`Lines` object.
The modifier is typically used to visualize the paths of motion of particles as continuous curves.
See the corresponding user manual page for more information.

```python
  from ovito.io import import_file
  from ovito.modifiers import GenerateTrajectoryLinesModifier
  
  # Load a particle simulation sequence:
  pipeline = import_file("input/simulation.*.dump")
  
  # Insert the modifier into the pipeline for creating the trajectory lines.
  # It will automatically sample the particle positions over the entire animation interval.
  modifier = GenerateTrajectoryLinesModifier(only_selected=False)
  pipeline.modifiers.append(modifier)
  
  # Configure trajectory line visualization:
  modifier.vis.width = 0.4
  modifier.vis.color = (1, 0, 0)
  modifier.vis.flat_shading = False
  
  # Insert pipeline into the scene to make the particles and
  # the trajectory lines visible in rendered images.
  pipeline.add_to_scene()
```"""
    frame_interval: Optional[Tuple[int, int]] = None
    'frame_interval() -> Optional[tuple[int, int]]\n\nThe animation frame interval over which the particle positions are sampled to generate the trajectory lines. Set this to a tuple of two integers to specify the first and the last animation frame; or use ``None`` to generate trajectory lines over the entire animation sequence.\n\nDefault: ``None``'
    only_selected: bool = True
    'only_selected() -> bool\n\nControls whether or not trajectory lines are generated only for selected particles. If this option is enabled, trajectory lines will be generated only for the subset of particles that are selected in the initial frame of the trajectory sequence, i.e., at the beginning of the :py:attr:`frame_interval`. If disabled, trajectory lines are generated for *all* particles. \n\nDefault: ``True``'
    sample_particle_property: str = ''
    'sample_particle_property() -> str\n\nName of the particle property to sample along the generated trajectory lines. \n\nThis option transfers the specified particle property to the trajectory line vertices. In other words, the time-dependent input property of the particles becomes a position-dependent property of the generated trajectory lines. The sampled values will become available in the :py:class:`Lines` object and may subsequently be used for pseudo-coloring of the trajectory lines. See the `LinesVis.color_mapping_property` option. \n\nThe following code example demonstrates how to color trajectory lines based on the velocity magnitude of the particles at the corresponding simulation times: \n\n\n```python\n  from ovito.io import *\n  from ovito.modifiers import *\n  from ovito.vis import *\n  import numpy\n  \n  # Import an MD simulation trajectory. The snapshots stored in the trajectory file contain the instantaneous\n  # particle velocities (vector property "Velocity") at each timestep and OVITO automatically\n  # computes the scalar velocities (particle property "Velocity Magnitude") during file import.\n  pipeline = import_file("../../../ovito/tests/files/GSD/trajectory.gsd")\n  \n  # Turn off the display of particles.\n  pipeline.compute().particles.vis.enabled = False\n  \n  # Configure modifier:\n  modifier = GenerateTrajectoryLinesModifier(only_selected=False)\n  modifier.sample_particle_property = "Velocity Magnitude"\n  # Configure trajectory line display:\n  modifier.vis.width = 0.1\n  modifier.vis.flat_shading = False\n  \n  # Locally color the trajectory lines based on the \'Velocity Magnitude\' property, which was transferred\n  # from the particles to the trajectory line vertices.\n  modifier.vis.color_mapping_property = "Velocity Magnitude"\n  \n  # Insert modifier into the pipeline to generate the trajectory lines.\n  pipeline.modifiers.append(modifier)\n  \n  # To determine the correct min/max value range for pseudo-coloring, we need to access the\n  # array of property values associated with the generated trajectory lines:\n  trajectories = pipeline.compute().lines["trajectories"]\n  \n  # Now we can set the min/max value range for the color mapping.\n  # Note that this change to the visual element does NOT trigger a re-evaluation of the modifier.\n  modifier.vis.color_mapping_interval = (\n      numpy.min(trajectories["Velocity Magnitude"]),\n      numpy.max(trajectories["Velocity Magnitude"]),\n  )\n  \n  # Insert pipeline into the scene to make the trajectory lines visible in rendered images.\n  pipeline.add_to_scene()\n```\n\n:emphasize-lines: 21-23\n\n\nDefault: ``\'\'``'
    sampling_frequency: int = 1
    'sampling_frequency() -> int\n\nLength of the animation frame intervals at which the particle positions should be sampled.\n\nDefault: ``1``'
    unwrap_trajectories: bool = True
    'unwrap_trajectories() -> bool\n\nControls whether or not trajectory lines should be automatically unwrapped at the box boundaries when the particles cross a periodic boundary.\n\nDefault: ``True``'
    vis: ovito.vis.LinesVis = ovito.vis.LinesVis()
    'vis() -> ovito.vis.LinesVis\n\nThe :py:class:`LinesVis` element controlling the visual appearance of the trajectory lines created by this modifier.'

@dataclass(kw_only=True)
class GrainSegmentationModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

This analysis modifier identifies the individual grains in a polycrystalline microstructure.
See the corresponding user manual page for 
more information on this modifier.

Inputs:

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Particle properties
      -
    * - ``Position``
      - The input coordinates of the particles.
    * - ``Structure Type``
      - The local structure types computed by the :py:class:`PolyhedralTemplateMatchingModifier`.
    * - ``Orientation``
      - The local lattice orientations computed by the :py:class:`PolyhedralTemplateMatchingModifier`.

Outputs:

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Particle properties
      -
    * - ``Grain``
      - The numeric identifier of the grain a particle has been assigned to by the algorithm.
    * - ``Color``
      - A per-particle color representing either the grain (if :py:attr:`color_particles` is ``True``) or the identified PTM structure type.

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Global attributes
      -
    * - ``GrainSegmentation.grain_count``
      - Number of grains that have been found by the algorithm.

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - :py:class:`DataTable`
      -
    * - ``grains``
      - A data table containing the list of grains found by the algorithm.
    * - ``grains-merge``
      - A scatter plot of the ordered cluster merge steps computed by the grain segmentation algorithm."""

    class Algorithm(enum.Enum):
        """"""
        GraphClusteringAuto = enum.auto()
        GraphClusteringManual = enum.auto()
        MinimumSpanningTree = enum.auto()
    algorithm: GrainSegmentationModifier.Algorithm = Algorithm.GraphClusteringAuto
    'algorithm() -> GrainSegmentationModifier.Algorithm\n\nSelects the cluster merge algorithm to use. One of the following constants:\n\n  * ``GrainSegmentationModifier.Algorithm.GraphClusteringAuto``\n  * ``GrainSegmentationModifier.Algorithm.GraphClusteringManual``\n  * ``GrainSegmentationModifier.Algorithm.MinimumSpanningTree``\n\n\nDefault: ``GrainSegmentationModifier.Algorithm.GraphClusteringAuto``'
    color_particles: bool = True
    'color_particles() -> bool\n\nControls whether the modifier colors atoms according to the grains they have been assigned to. \n\nDefault: ``True``'
    handle_stacking_faults: bool = True
    'handle_stacking_faults() -> bool\n\nControls whether the algorithm should handle stacking faults and coherent twin boundaries explicitly. If turned off, atoms with cubic and hexagonal crystal structure will simply be treated as separate phases and will never be merged into the same grain. \n\nDefault: ``True``'
    merging_threshold: float = 0.0
    'merging_threshold() -> float\n\nSpecifies the maximum graph edge contraction distance and determines the resulting number and sizes of grains. A lower threshold produces more (and smaller) grains; a larger threshold produces fewer (and larger) grains. \n\nThis parameter is ignored if :py:attr:`algorithm` is ``GraphClusteringAuto``, in which case the merging threshold is automatically determined by the algorithm to give optimal segmentation results. \n\nDefault: ``0.0``'
    min_grain_size: int = 100
    'min_grain_size() -> int\n\nMinimum number of atoms in a grain. Crystallite clusters smaller than this threshold get dissolved during the grain segmentation and their atoms are merged into neighboring grains. \n\nDefault: ``100``'
    orphan_adoption: bool = True
    'orphan_adoption() -> bool\n\nControls the merging of non-crystalline atoms (e.g. grain boundary atoms) into adjacent grains. \n\nDefault: ``True``'

@dataclass(kw_only=True)
class HistogramModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

Generates a histogram of a property, i.e. the distribution of its per-element values.
See also the corresponding user manual page for more information.
The modifier can operate on different types of data elements:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Data element
      -
    * - ``particles``
      - Computes the histogram for a particle property.
    * - ``bonds``
      - Computes the histogram for a bond property.
    * - ``voxels``
      - Computes the histogram for a voxel grid property.

By default the modifier will operate on a particle property. You can change this by setting the :py:attr:`operate_on` field.

The value range of the histogram is determined automatically from the minimum and maximum values of the selected property
unless :py:attr:`fix_xrange` is set to ``True``. In this case the range of the histogram is controlled by the
:py:attr:`xrange_start` and :py:attr:`xrange_end` parameters.

Example:

```python
  from ovito.modifiers import HistogramModifier
  from ovito.io import import_file, export_file
  
  pipeline = import_file("input/simulation.dump")
  
  modifier = HistogramModifier(bin_count=100, property='peatom')
  pipeline.modifiers.append(modifier)
  
  export_file(pipeline, "output/histogram.txt", "txt/table", key="histogram[peatom]")
```"""
    bin_count: int = 200
    'bin_count() -> int\n\nThe number of histogram bins.\n\nDefault: ``200``'
    fix_xrange: bool = False
    "fix_xrange() -> bool\n\nControls how the value range of the histogram is determined. If set to ``False``, the interval and bin widths are chosen automatically by the modifier to include all input values. If set to ``True``, the histogram's value range may be specified explicitly by setting :py:attr:`xrange_start` and :py:attr:`xrange_end`.\n\nDefault: ``False``"
    only_selected: bool = False
    'only_selected() -> bool\n\nIf ``True``, the histogram is computed only on the basis of currently selected particles or bonds. You can use this to restrict histogram calculation to a subset of input particles/bonds. \n\nDefault: ``False``'
    operate_on: str = 'particles'
    "operate_on() -> str\n\nSelects the kind of data elements this modifier should operate on. Supported values are: ``'particles'``, ``'bonds'``, ``'voxels'``. \n\nDefault: ``'particles'``"
    property: str = ''
    'property() -> str\n\nThe name of the input property for which to compute the histogram. For vector properties a component name must be appended in the string, e.g. ``"Velocity.X"``. \n\nDefault: ``\'\'``'
    xrange_end: float = 0.0
    'xrange_end() -> float\n\nIf :py:attr:`fix_xrange` is true, then this parameter controls the upper end of the value interval covered by the histogram. Input values higher than this range limit will be discarded from the computed histogram. \n\nDefault: ``0.0``'
    xrange_start: float = 0.0
    'xrange_start() -> float\n\nIf :py:attr:`fix_xrange` is true, then this parameter controls the lower end of the value interval covered by the histogram.Input values lower than this range limit will be discarded from the computed histogram. \n\nDefault: ``0.0``'
    select_elements: bool = False
    'select_elements() -> bool\n\nLets the modifier create a selection containing all elements whose property values are in the range [:py:attr:`selection_start`, :py:attr:`selection_end`]. \n\nDefault: ``False``\n\n\n.. versionadded: 3.11.0'
    selection_start: float = 0.0
    'selection_start() -> float\n\nLower bound of the selection interval if the :py:attr:`select_elements` option is turned on.\n\nDefault: ``0.0``\n\n\n.. versionadded: 3.11.0'
    selection_end: float = 0.0
    'selection_end() -> float\n\nUpper bound of the selection interval if the :py:attr:`select_elements` option is turned on.\n\nDefault: ``0.0``\n\n\n.. versionadded: 3.11.0'

@dataclass(kw_only=True)
class IdentifyDiamondModifier(StructureIdentificationModifier):
    """Base: :py:class:`ovito.modifiers.StructureIdentificationModifier`

Analyzes the local neighborhood of each particle to identify cubic or hexagonal diamond lattices.
See the corresponding user manual page
for more information.

Note that this class inherits several important parameter fields from its :py:class:`StructureIdentificationModifier`
base class.

Modifier inputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Position``
      - The input coordinates of the particles.
    * - ``Selection``
      - The selection state of the input particles. Only needed if :py:attr:`~StructureIdentificationModifier.only_selected` is ``True``.

Modifier outputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Structure Type``
      - The structure type computed by the algorithm for each particle, encoded as an integer value:

        ============= =========================================================
        Numeric id    Python constant
        ============= =========================================================
        0             ``IdentifyDiamondModifier.Type.OTHER``
        1             ``IdentifyDiamondModifier.Type.CUBIC_DIAMOND``
        2             ``IdentifyDiamondModifier.Type.CUBIC_DIAMOND_FIRST_NEIGHBOR``
        3             ``IdentifyDiamondModifier.Type.CUBIC_DIAMOND_SECOND_NEIGHBOR``
        4             ``IdentifyDiamondModifier.Type.HEX_DIAMOND``
        5             ``IdentifyDiamondModifier.Type.HEX_DIAMOND_FIRST_NEIGHBOR``
        6             ``IdentifyDiamondModifier.Type.HEX_DIAMOND_SECOND_NEIGHBOR``
        ============= =========================================================
    * - ``Color``
      - Particle coloring to indicate the identified structure type for each particle; only if :py:attr:`~StructureIdentificationModifier.color_by_type` is ``True``.
        See the :py:attr:`~StructureIdentificationModifier.structures` array on how to customize the colors.

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Global attributes
      -
    * - ``IdentifyDiamond.counts.OTHER``
      - Number of particles not matching any of the recognized structure types.
    * - ``IdentifyDiamond.counts.CUBIC_DIAMOND``
      - Number of particles identified as fully coordinated cubic diamond.
    * - ``IdentifyDiamond.counts.CUBIC_DIAMOND_FIRST_NEIGHBOR``
      - Number of particles identified as partially coordinated cubic diamond.
    * - ``IdentifyDiamond.counts.CUBIC_DIAMOND_SECOND_NEIGHBOR``
      - Number of particles identified as partially coordinated cubic diamond.
    * - ``IdentifyDiamond.counts.HEX_DIAMOND``
      - Number of particles identified as fully coordinated hexagonal diamond.
    * - ``IdentifyDiamond.counts.HEX_DIAMOND_FIRST_NEIGHBOR``
      - Number of particles identified as partially coordinated hexagonal diamond.
    * - ``IdentifyDiamond.counts.HEX_DIAMOND_SECOND_NEIGHBOR``
      - Number of particles identified as partially coordinated hexagonal diamond.

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Data tables
      -
    * - ``structures``
      - A bar chart with the particle counts for each structure type identified by the modifier.
        You can retrieve this :py:class:`DataTable` from the `DataCollection.tables` dictionary."""

    class Type(enum.IntEnum):
        """"""
        OTHER = enum.auto()
        CUBIC_DIAMOND = enum.auto()
        CUBIC_DIAMOND_FIRST_NEIGHBOR = enum.auto()
        CUBIC_DIAMOND_SECOND_NEIGHBOR = enum.auto()
        HEX_DIAMOND = enum.auto()
        HEX_DIAMOND_FIRST_NEIGHBOR = enum.auto()
        HEX_DIAMOND_SECOND_NEIGHBOR = enum.auto()

@dataclass(kw_only=True)
class IdentifyFCCPlanarFaultsModifier(ovito.pipeline.ModifierInterface):
    """Base: :py:class:`ovito.pipeline.ModifierInterface`

This Python-based modifier identifies different planar defect types in face-centered cubic (fcc) crystals
such as stacking faults and coherent twin boundaries. All of these planar defects have in common that they are
made of hcp-like atoms arranged in {111} planes of the fcc crystal.
See also the corresponding user manual page for more information.

The identification algorithm relies on intermediate results provided by the :py:class:`PolyhedralTemplateMatchingModifier`, which
must be inserted into the data pipeline first to identify all hcp-like defect atoms in the fcc crystal:

```python
  from ovito.modifiers import PolyhedralTemplateMatchingModifier, IdentifyFCCPlanarFaultsModifier
  pipeline.modifiers.append(PolyhedralTemplateMatchingModifier(output_orientation=True, output_interatomic_distance=True))
  pipeline.modifiers.append(IdentifyFCCPlanarFaultsModifier())
```

The modifier subsequently analyzes the neighborhood of each hcp-like atom to identify which kind
of planar fault it is part of. Each atom in the input system is attributed to one of following classes:

  * 0 = Non-hcp structure (e.g. perfect fcc)
  * 1 = Indeterminate hcp-like (e.g. isolated hcp-like atoms, not a planar defect)
  * 2 = Intrinsic stacking fault (two hcp-like layers)
  * 3 = Coherent twin boundary (one hcp-like layer)
  * 4 = Multilayer stacking fault (three or more hcp-like layers)

The modifier writes the corresponding numeric values to the ``Planar Fault Type`` output particle property.
Furthermore, the modifier produces a :py:class:`DataTable` which lists the total atom count and
estimated total area for each planar defect type. This information can be accessed as follows:

```python
  data = pipeline.compute()
  table = data.tables['planar_faults']
  num_atoms_isf  = table['Atom Count'][IdentifyFCCPlanarFaultsModifier.Types.ISF]
  num_atoms_twin = table['Atom Count'][IdentifyFCCPlanarFaultsModifier.Types.TWIN]
  area_isf       = table['Estimated Area'][IdentifyFCCPlanarFaultsModifier.Types.ISF]
  area_twin      = table['Estimated Area'][IdentifyFCCPlanarFaultsModifier.Types.TWIN]
```"""

    class Types(enum.IntEnum):
        """The possible classifications for hcp-like defect atoms."""
        NONHCP = 0
        OTHER = 1
        ISF = 2
        TWIN = 3
        MULTI = 4
    color_other: ovito.vis.Color = (1.0, 0.4, 0.4)
    color_isf: ovito.vis.Color = (0.9, 0.7, 0.2)
    color_twin: ovito.vis.Color = (0.4, 0.4, 1.0)
    color_multi: ovito.vis.Color = (0.0, 1.0, 0.9)

@dataclass(kw_only=True)
class InvertSelectionModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

This modifier inverts the current data element selection.
See also the corresponding user manual page for more information.
The modifier can operate on different kinds of data elements:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Data element
      -
    * - ``particles``
      - Inverts the values of the ``Selection`` particle property.
    * - ``bonds``
      - Inverts the values of the ``Selection`` bond property.
    * - ``voxels``
      - Inverts the values of the ``Selection`` voxel grid property.

By default the modifier will act on particles. You can change this by setting the :py:attr:`operate_on` field."""
    operate_on: str = 'particles'
    "operate_on() -> str\n\nSelects the kind of data elements this modifier should operate on. Supported values are: ``'particles'``, ``'bonds'``, ``'voxels'``. \n\nDefault: ``'particles'``"

@dataclass(kw_only=True)
class LoadTrajectoryModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

This modifier loads trajectories of particles from a separate simulation file.
See also the corresponding user manual page for this modifier.

The typical use case for this modifier is when the topology of a molecular system (i.e. the definition of atom types, bonds, etc.) is
stored separately from the trajectories of atoms. In this case you should load the topology file first using :py:func:`import_file`.
Then create and apply the :py:class:`LoadTrajectoryModifier` to the topology dataset, which loads the trajectory file.
The modifier will replace the static atom positions from the topology dataset with the time-dependent positions from the trajectory file.

Example:

```python
  from ovito.io import import_file
  from ovito.modifiers import LoadTrajectoryModifier
  
  # Load static topology data from a LAMMPS data file.
  pipeline = import_file('input/input.data', atom_style='bond')
  
  # Load atom trajectories from separate LAMMPS dump file.
  traj_mod = LoadTrajectoryModifier()
  traj_mod.source.load('input/trajectory.dump')
  print("Number of frames: ", traj_mod.source.num_frames)
  
  # Insert modifier into data pipeline.
  pipeline.modifiers.append(traj_mod)
```

Furthermore, it is possible to let the modifier load varying bond topologies
for each trajectory frame from a LAMMPS `dump local <https://docs.lammps.org/dump.html>`__ file:

```python
  bonds_mod = LoadTrajectoryModifier()
  bonds_mod.source.load('input/bonds.dump.local', 
      columns = [None, 'Bond Type', 'Particle Identifiers.1', 'Particle Identifiers.2', 'Length', 'Energy'])
  pipeline.modifiers.append(bonds_mod)
```

Here, the ``columns`` function parameter specifies the mapping of data columns in the imported dump file to corresponding
target bond properties within OVITO. The dump local file :file:`bonds.dump.local` contains six data columns
and has been produced by the following LAMMPS commands::

   compute 1 all property/local btype batom1 batom2
   compute 2 all bond/local dist engpot
   dump bonds all local 100 bonds.dump.local index c_1[1] c_1[2] c_1[3] c_2[1] c_2[2]"""
    source: ovito.pipeline.FileSource = ovito.pipeline.FileSource()
    'source() -> ovito.pipeline.FileSource\n\nA :py:class:`FileSource` that provides the trajectories of particles. You can call its :py:meth:`load` function to load a simulation trajectory file as shown in the code example above.'

@dataclass(kw_only=True)
class PolyhedralTemplateMatchingModifier(StructureIdentificationModifier):
    """Base: :py:class:`ovito.modifiers.StructureIdentificationModifier`

This modifier analyzes the local neighborhood of each particle to identify common structural motives and crystalline structures.
The structure identification is based on the Polyhedral Template Matching (PTM) algorithm.
See the corresponding user manual page
for more information. The PTM algorithm is able to compute local crystal orientations, elastic lattice strains, and can identify 
local chemical orderings in binary compounds.

Note that this modifier inherits several important parameter fields from the :py:class:`StructureIdentificationModifier`
base class.

Modifier inputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Position``
      - The input coordinates of the particles.
    * - ``Particle Type``
      - The chemical types of the input particles; only used if :py:attr:`output_ordering` is ``True``.
    * - ``Selection``
      - The selection state of the input particles; only used if :py:attr:`~StructureIdentificationModifier.only_selected` is ``True``.

Modifier outputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Structure Type``
      - The structure type computed by the algorithm for each particle, encoded as an integer value:

        ============= ========================================================= =============
        Numeric id    Python constant                                           Initial state
        ============= ========================================================= =============
        0             ``PolyhedralTemplateMatchingModifier.Type.OTHER``         
        1             ``PolyhedralTemplateMatchingModifier.Type.FCC``           enabled
        2             ``PolyhedralTemplateMatchingModifier.Type.HCP``           enabled
        3             ``PolyhedralTemplateMatchingModifier.Type.BCC``           enabled
        4             ``PolyhedralTemplateMatchingModifier.Type.ICO``           disabled
        5             ``PolyhedralTemplateMatchingModifier.Type.SC``            disabled
        6             ``PolyhedralTemplateMatchingModifier.Type.CUBIC_DIAMOND`` disabled
        7             ``PolyhedralTemplateMatchingModifier.Type.HEX_DIAMOND``   disabled
        8             ``PolyhedralTemplateMatchingModifier.Type.GRAPHENE``      disabled
        ============= ========================================================= =============

        The algorithm only identifies enabled structure types; see :py:attr:`~StructureIdentificationModifier.structures` array for details.
    * - ``RMSD``
      - The per-particle RMSD values computed by the PTM algorithm.
        Only if :py:attr:`output_rmsd` is set.
    * - ``Interatomic Distance``
      - The per-particle local atomic distances computed by the PTM algorithm.
        Only if :py:attr:`output_interatomic_distance` is set.
    * - ``Orientation``
      - The local lattice orientations computed by the PTM algorithm, encoded as quaternions.
        Only if :py:attr:`output_orientation` is set.
    * - ``Elastic Deformation Gradient``
      - The per-particle elastic deformation gradient tensors computed by the PTM algorithm (3x3 components).
        Only if :py:attr:`output_deformation_gradient` is set.
    * - ``Ordering Type``
      - The local chemical ordering type determined by the PTM algorithm, encoded as an integer value.
        Only if :py:attr:`output_ordering` is set.

        ============= =========================================================
        Numeric id    Python constant
        ============= =========================================================
        0             ``PolyhedralTemplateMatchingModifier.OrderingType.OTHER``
        1             ``PolyhedralTemplateMatchingModifier.OrderingType.PURE``
        2             ``PolyhedralTemplateMatchingModifier.OrderingType.L10``
        3             ``PolyhedralTemplateMatchingModifier.OrderingType.L12_A``
        4             ``PolyhedralTemplateMatchingModifier.OrderingType.L12_B``
        5             ``PolyhedralTemplateMatchingModifier.OrderingType.B2``
        6             ``PolyhedralTemplateMatchingModifier.OrderingType.ZINCBLENDE_WURTZITE``
        7             ``PolyhedralTemplateMatchingModifier.OrderingType.BORON_NITRIDE``
        ============= =========================================================
    * - ``Color``
      - Particle coloring to indicate the identified structure type for each particle; only if :py:attr:`~StructureIdentificationModifier.color_by_type` is ``True``.
        See the :py:attr:`~StructureIdentificationModifier.structures` array on how to customize the colors.

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Global attributes
      -
    * - ``PolyhedralTemplateMatching.counts.OTHER``
      - Number of particles not matching any of the recognized structural types.
    * - ``PolyhedralTemplateMatching.counts.FCC``
      - Number of particles identified as face-centered cubic structure.
    * - ``PolyhedralTemplateMatching.counts.HCP``
      - Number of particles identified as hexagonal close packed structure.
    * - ``PolyhedralTemplateMatching.counts.BCC``
      - Number of particles identified as body-centered cubic structure.
    * - ``PolyhedralTemplateMatching.counts.ICO``
      - Number of particles identified as icosahedral structure.
    * - ``PolyhedralTemplateMatching.counts.SC``
      - Number of particles identified as simple cubic structure.
    * - ``PolyhedralTemplateMatching.counts.CUBIC_DIAMOND``
      - Number of particles identified as cubic diamond structure.
    * - ``PolyhedralTemplateMatching.counts.HEX_DIAMOND``
      - Number of particles identified as hexagonal diamond structure.
    * - ``PolyhedralTemplateMatching.counts.GRAPHENE``
      - Number of particles identified as 2d graphene structure.

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Data tables
      -
    * - ``structures``
      - A bar chart with the particle counts for each structure type identified by the modifier.
        You can retrieve this :py:class:`DataTable` from the `DataCollection.tables` dictionary.
    * - ``ptm-rmsd``
      - A histogram with the RMSD value distribution computed by the modifier.
        You can retrieve this :py:class:`DataTable` from the `DataCollection.tables` dictionary."""

    class Type(enum.IntEnum):
        """"""
        OTHER = enum.auto()
        FCC = enum.auto()
        HCP = enum.auto()
        BCC = enum.auto()
        ICO = enum.auto()
        SC = enum.auto()
        CUBIC_DIAMOND = enum.auto()
        HEX_DIAMOND = enum.auto()
        GRAPHENE = enum.auto()

    class OrderingType(enum.IntEnum):
        """"""
        NONE = enum.auto()
        PURE = enum.auto()
        L10 = enum.auto()
        L12_A = enum.auto()
        L12_B = enum.auto()
        B2 = enum.auto()
        ZINCBLENDE_WURTZITE = enum.auto()
        BORON_NITRIDE = enum.auto()
    output_deformation_gradient: bool = False
    'output_deformation_gradient() -> bool\n\nBoolean flag that controls whether the modifier outputs the computed per-particle elastic deformation gradients as a new particle property named ``Elastic Deformation Gradient``.The elastic deformation gradient describes the local deformation and rigid-body rotation of the crystal with respect to an ideal reference lattice configuration. See the OVITO user manual for details. \n\nDefault: ``False``'
    output_interatomic_distance: bool = False
    'output_interatomic_distance() -> bool\n\nBoolean flag that controls whether the modifier outputs the computed per-particle interatomic distance as a new particle property named ``Interatomic Distance``.\n\nDefault: ``False``'
    output_ordering: bool = False
    'output_ordering() -> bool\n\nBoolean flag that controls whether the modifier should identify local chemical ordering types and output them as a new particle property named ``Ordering Type``. \n\nDefault: ``False``'
    output_orientation: bool = False
    'output_orientation() -> bool\n\nBoolean flag that controls whether the modifier outputs the computed per-particle lattice orientations as a new particle property named ``Orientation``. The lattice orientation is specified in terms of a quaternion that describes the rotation of the crystal with respect to a reference lattice orientation. See the OVITO user manual for details. \n\nDefault: ``False``'
    output_rmsd: bool = False
    'output_rmsd() -> bool\n\nBoolean flag that controls whether the modifier outputs the computed per-particle RMSD values as a new particle property named ``RMSD``.\n\nDefault: ``False``'
    rmsd_cutoff: float = 0.1
    'rmsd_cutoff() -> float\n\nThe maximum allowed root mean square deviation for positive structure matches. If this threshold value is non-zero, template matches that yield a RMSD value above the cutoff are classified as "Other". This can be used to filter out spurious template matches (false positives). \n\nDefault: ``0.1``'

    @staticmethod
    def calculate_misorientation(a: ArrayLike, b: ArrayLike, symmetry: str='none', output: str='angle') -> NDArray[Any]:
        """Helper function for computing the `misorientation or disorientation <https://en.wikipedia.org/wiki/Misorientation>`__ of a pair of crystal orientations in 3d space. The *misorientation* is defined as the 3d rotation mapping crystal orientation *a* to orientation *b* -- without taking into account lattice symmetries. The *disorientation* is obtained by subsequently mapping the misorientation into the fundamental zone, i.e., minimizing the rotation angle by considering all equivalent orientations due to lattice symmetry. 

Orientations *a* and *b* must both be specified as quaternions of the form :math:`(x,y,z,w)`, i.e., as arrays of size 4. :py:class:`PolyhedralTemplateMatchingModifier` outputs the local lattice orientations of crystalline atoms in this form (see :py:attr:`output_orientation`) and :py:class:`GrainSegmentationModifier` calculates the mean lattice orientation of each grain in this form for a polycrystal. 

The parameter *output* controls whether the function should return just the scalar angle of rotation (in radians) or the full 3d crystal rotation (as a quaternion). 

:param a: First crystal orientation in the form of a 4-component quaternion.
:param b: Second crystal orientation in the form of a 4-component quaternion.
:param symmetry: Selects between *misorientation* and *disorientation* calculation:

   * ``'none'``: Calculate misorientation (ignore lattice symmetry).
   * ``'cubic'``: Calculate disorientation for lattices with cubic symmetry (e.g. fcc, bcc, diamond).
   * ``'hexagonal'``: Calculate disorientation for lattices with hexagonal symmetry (e.g. hcp, hex-diamond, graphene).
:param output: Controls the type of output value returned by the function:

   * ``'angle'``: Return scalar misorientation/disorientation angle in radians.
   * ``'rotation'``: Return 3d rotation from *a* to *b* expressed as quaternion.


Parameters *a* and *b* may both be arrays of shape (*N*,4) specifying two orientation lists of the same length *N*. In this case the function will compute the list of *N* pair-wise misorientations in one go, which is usually much faster than calling the function in a loop. 

Usage example:

```python
  pipeline.modifiers.append(PolyhedralTemplateMatchingModifier(output_orientation=True))
  pipeline.modifiers.append(GrainSegmentationModifier())
  data = pipeline.compute()
  
  # Assuming the input structure is an fcc bicrystal consisting of exactly two grains.
  assert data.attributes['GrainSegmentation.grain_count'] == 2 
  grain_orientations = data.tables['grains']['Orientation']
  
  # Calculate the disorientation angle between the two bicrystal grains.
  disorientation_angle = PolyhedralTemplateMatchingModifier.calculate_misorientation(
      grain_orientations[0], grain_orientations[1], 
      symmetry='cubic')
```"""
        ...

@dataclass(kw_only=True)
class PythonModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

This is a wrapper for user-defined Python modifiers, which allows inserting them into a :py:class:`Pipeline`.

This class makes it possible to implement new modifier functions in the Python language which can participate in OVITO's
data pipeline system and which can be used just like OVITO's built-in modifiers.
Please see the :ref:`writing_custom_modifiers` section.

Example using the simple programming interface:

```python
  from ovito.io import import_file
  
  # Load some input data:
  pipeline = import_file("input/simulation.dump")
  
  # Define our custom modifier function, which assigns a uniform color
  # to all particles, similar to the built-in AssignColorModifier.
  def assign_color(frame, data):
      color_property = data.particles_.create_property('Color')
      color_property[:] = (0.2, 0.5, 1.0)
  
  # Insert the user-defined modifier function into the data pipeline.
  pipeline.modifiers.append(assign_color)
  # Note that appending the Python function to the pipeline is equivalent to
  # creating a PythonModifier instance and appending it:
  #
  #   pipeline.modifiers.append(PythonModifier(function=assign_color))
  
  # Evaluate data pipeline. This will make the system invoke assign_color().
  data = pipeline.compute()
  print(data.particles.colors[...])
```

Example using the advanced programming interface:

```python
  from ovito.io import import_file
  from ovito.pipeline import ModifierInterface
  
  # Load some input data:
  pipeline = import_file("input/simulation.dump")
  
  # Define our custom modifier class, which assigns a uniform color
  # to all particles, similar to the built-in AssignColorModifier.
  class AssignColor(ModifierInterface):
  
      def modify(self, data, **kwargs):
          color_property = data.particles_.create_property('Color')
          color_property[:] = (0.2, 0.5, 1.0)
  
  # Insert the user-defined modifier into the data pipeline.
  pipeline.modifiers.append(AssignColor())
  # Note that appending an instance of your class to the pipeline is equivalent to
  # wrapping it in a PythonModifier instance:
  #
  #   pipeline.modifiers.append(PythonModifier(delegate=AssignColor()))
  
  # Evaluate data pipeline. This will make the system invoke AssignColor.modify().
  data = pipeline.compute()
  print(data.particles.colors[...])
```"""
    function: Optional[Callable[[int, ovito.data.DataCollection], Optional[Generator[str | float, None, None]]]] = None
    'function() -> Optional[collections.abc.Callable[[int, ovito.data.DataCollection], Any]]\n\nThe Python function to be called each time the data pipeline is evaluated by the system. This field is only used if the user-defined modifier is implemented using the simple programming interface. \n\nThe function must have a particular signature as shown in the first example above. The *frame* parameter specifies the current animation frame number at which the data pipeline is being evaluated. The :py:class:`DataCollection` *data* initially holds the input data objects of the modifier, which were produced by the upstream part of the data pipeline. The user-defined modifier function is free modify the data collection and the data objects stored in it. \n\nDefault: ``None``'
    delegate: Optional[ovito.pipeline.ModifierInterface] = None
    'delegate() -> Optional[ovito.pipeline.ModifierInterface]\n\nA :py:class:`ModifierInterface` object implementing the logic of the user-defined modifier. This field is only used if the user-defined modifier is implemented using the advanced programming interface. \n\nDefault: ``None``'
    working_dir: str = ''
    "working_dir() -> str\n\nA path that will be set as active working directory while the Python modifier function is executed by the pipeline system. This setting mainly plays a role if the modifier function is used within the GUI of OVITO and if it performs some file I/O. Relative file paths will then get resolved with respect to this working directory. \n\nIf no specific working directory is set for the modifier function, the application's current working directory will be used. \n\nDefault: ``''``"

@dataclass(kw_only=True)
class RenderLAMMPSRegionsModifier(ovito.pipeline.ModifierInterface):
    """Base: :py:class:`ovito.pipeline.ModifierInterface`

This Python-based modifier allows visualizing
spatial regions with different 3d geometries as defined by the `region <https://docs.lammps.org/region.html>`__ command
of the `LAMMPS <https://docs.lammps.org/>`__ simulation code.
See also the corresponding user manual page for more information on this OVITO modifier.

The :py:class:`SurfaceMesh` generated by the modifier is accessible
in the pipeline output :py:class:`DataCollection` under the unique :py:attr:`identifier` *lammps-regions*::

    data = pipeline.compute()
    mesh = data.surfaces['lammps-regions']"""
    commands: str = 'region 1 sphere -20 0 0 15\nregion 2 cylinder z 10 0 5 INF INF'
    select_atoms: bool = False
    color_by_region: bool = False
    vis: ovito.vis.SurfaceMeshVis = ovito.vis.SurfaceMeshVis(smooth_shading=False, clip_at_domain_boundaries=True, title='LAMMPS Regions')

@dataclass(kw_only=True)
class ReplicateModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

This modifier replicates all particles, bonds and other elements of a system to visualize periodic images.
See also the corresponding user manual page for more information.

Inputs:

The modifier can operate on any combination of the following data elements:

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Data element specifier
      - Description
    * - ``particles``
      - Duplicates :py:class:`Particles` and :py:class:`Bonds`.
    * - ``surfaces``
      - Duplicates the mesh geometry of :py:class:`SurfaceMesh` objects.
    * - ``voxels``
      - Duplicates the voxel elements of :py:class:`VoxelGrid` objects.
    * - ``vectors``
      - Transforms the base points and directions of all :py:class:`Vectors` objects.
    * - ``lines``
      - Duplicates the lines stored in all :py:class:`Lines` objects.
    * - ``dislocations``
      - Duplicates dislocation lines in a :py:class:`DislocationNetwork`.

By default the modifier will operate on all of these. You can restrict it to a subset by setting the :py:attr:`operate_on` field."""
    adjust_box: bool = True
    'adjust_box() -> bool\n\nControls whether the simulation cell is resized. If ``True``, the simulation cell is accordingly extended to fit the replicated data. If ``False``, the original simulation cell size (containing only the primary image of the system) is maintained. \n\nDefault: ``True``'
    num_x: int = 1
    'num_x() -> int\n\nA positive integer specifying the number of copies to generate in the *x* direction (including the existing primary image).\n\nDefault: ``1``'
    num_y: int = 1
    'num_y() -> int\n\nA positive integer specifying the number of copies to generate in the *y* direction (including the existing primary image).\n\nDefault: ``1``'
    num_z: int = 1
    'num_z() -> int\n\nA positive integer specifying the number of copies to generate in the *z* direction (including the existing primary image).\n\nDefault: ``1``'
    operate_on: MutableSet[str] = {'particles', 'voxels', 'surfaces', 'lines', 'dislocations'}
    "A set of strings specifying the kinds of data elements this modifier should operate on. By default the set contains all data element types supported by the modifier. \n\nDefault: ``{'particles', 'dislocations', 'lines', 'surfaces','vectors', 'voxels'}``"
    unique_ids: bool = True
    "unique_ids() -> bool\n\nIf ``True``, the modifier automatically generates new unique IDs for each copy of particles. Otherwise, the replica will keep the same IDs as the original particles, which is typically not what you want. \n\nNote: This option has no effect if the input particles do not already have numeric IDs (i.e. the ``'Particle Identifier'`` property does not exist). \n\nDefault: ``True``"

@dataclass(kw_only=True)
class SelectTypeModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

Selects data elements of a certain type or types (e.g. all atoms of a chemical species).
See also the corresponding user manual page for more information.
The modifier can operate on different data elements:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Data element
      -
    * - ``particles``
      - Selects all particles of a certain type.
    * - ``bonds``
      - Selects all bonds of a certain type.

By default the modifier will act on particles. You can change this by setting the :py:attr:`operate_on` field.

Outputs:

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Global attributes
      -
    * - ``SelectType.num_selected``
      - The number of data elements (particles/bonds) that have been selected by the modifier.

Example:

```python
  from ovito.io import import_file
  from ovito.modifiers import SelectTypeModifier, CommonNeighborAnalysisModifier
  
  pipeline = import_file("input/simulation.dump")
  
  # Insert a CNA modifier to determine the structural type of each atom:
  pipeline.modifiers.append(CommonNeighborAnalysisModifier())
  
  # Apply the SelectTypeModifier to select all atoms of FCC and HCP type:
  pipeline.modifiers.append(SelectTypeModifier(
      operate_on = "particles",
      property = "Structure Type",
      types = { CommonNeighborAnalysisModifier.Type.FCC,
                CommonNeighborAnalysisModifier.Type.HCP }
  ))
  
  # The SelectTypeModifier reports the number of selected elements as an attribute:
  data = pipeline.compute()
  print("Number of FCC/HCP atoms: %i" % data.attributes['SelectType.num_selected'])
```"""
    operate_on: str = 'particles'
    "operate_on() -> str\n\nControls the kind of data elements this modifier should select. Supported values are: ``'particles'``, ``'bonds'``. \n\nDefault: ``'particles'``"
    property: str = 'Particle Type'
    "property() -> str\n\nThe name of the property to use as input; must be an integer property. \n\nFor selecting particles, possible input properties are ``'Particle Type'`` and ``'Structure Type'``, for example. For selecting bonds, ``'Bond Type'`` is a typical input property. \n\nDefault: ``'Particle Type'``"
    types: AbstractSet[Union[str, int]] = set([])
    "Specifies the types to select. You can add numeric type *IDs* or type *names* to this set.\nType names will automatically be translated into corresponding numeric type IDs by the modifier.\nThus, it is not necessary for you to look up the numeric ID for a type name using `Property.type_by_name()`.\nFor example, to select all atoms belonging to the type named 'Cu':\n\n```python\n  modifier = SelectTypeModifier(property = 'Particle Type', types = {'Cu'})\n```\n\nWhen using the :py:class:`SelectTypeModifier` to select *structural* types, you can directly use the predefined numeric constants for the structures, e.g.:\n\n```python\n  # Let the CNA modifier identify the structural type of each particle:\n  pipeline.modifiers.append(CommonNeighborAnalysisModifier())\n  # Select all FCC and HCP particles:\n  modifier = SelectTypeModifier(property = 'Structure Type')\n  modifier.types = { CommonNeighborAnalysisModifier.Type.FCC,\n                     CommonNeighborAnalysisModifier.Type.HCP }\n  pipeline.modifiers.append(modifier)\n```\n\nDefault: ``set([])``"

@dataclass(kw_only=True)
class SliceModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

Deletes or selects data elements located within a semi-infinite region bounded by a plane or within a slab bounded by a pair of parallel planes.
See also the corresponding user manual page for more information.

Inputs:

The modifier can operate on any combination of the following data elements:

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Data element
      -
    * - ``particles``
      - Deletes (or selects) particles.
    * - ``surfaces``
      - Cuts the mesh geometry of a :py:class:`SurfaceMesh`.
    * - ``voxels``
      - Creates a cross-sectional cut of a :py:class:`VoxelGrid`.
    * - ``lines``
      - Cuts the lines stored in :py:class:`Lines` objects.
    * - ``dislocations``
      - Cuts dislocation lines of a :py:class:`DislocationNetwork`.
    * - ``vectors``
      - Transforms the base points and directions of all :py:class:`Vectors` objects.

By default the modifier will act on all of these. You can restrict it to a subset by setting the :py:attr:`operate_on` field.
Furthermore, you can restrict the operation to only selected particles by setting the :py:attr:`only_selected` option.

The :py:attr:`select` option lets you select all elements on one side of the plane instead of deleting them.
Currently, the selection mode works only for :py:class:`Particles` or :py:class:`SurfaceMesh` vertices,
which get selected by setting their ``Selection`` particle property to 1."""
    distance: float = 0.0
    'distance() -> float\n\nThe distance of the slicing plane from the origin (along its normal vector).Its interpretation depends on whether :py:attr:`miller` mode is enabled or not. \n\nDefault: ``0.0``'
    gridslice_vis: ovito.vis.SurfaceMeshVis = ovito.vis.SurfaceMeshVis(show_cap=False, smooth_shading=False, title='Volume slice')
    "The :py:class:`SurfaceMeshVis` controlling the visual appearance of the cross-sectional slice being extracted from a :py:class:`VoxelGrid` by the modifier. The visual element is only used if :py:attr:`operate_on` includes ``'voxels'`` and the input data collection contains a :py:class:`VoxelGrid`."
    inverse: bool = False
    'inverse() -> bool\n\nReverses the sense of the slicing plane.\n\nDefault: ``False``'
    miller: bool = False
    'miller() -> bool\n\nControls whether the :py:attr:`.normal` vector and the :py:attr:`.distance` parameter are given in terms of the reciprocal cell vectors. \n\nIf enabled, the modifier will interpret the components of the :py:attr:`.normal` vector as Miller indices :math:`hkl`. Note that the indices do not have to be integer. The :attr:`.distance` parameter measures the (signed) offset of the plane from the simulation cell origin. It is specified in terms of the interplanar spacing :math:`d_{\\mathrm{hkl}}`, which depends on the simulation cell vectors and the Miller indices :math:`hkl`. \n\nIf Miller index mode is off, the :py:attr:`.normal` vector is specified in the Cartesian simulation coordinate system. The :attr:`.distance` from the origin (0,0,0) is measured in simulation units of length along the normal vector. \n\nThe :py:attr:`.slab_width` parameter is always specified in real-space units. \n\nDefault: ``False``'
    normal: ovito.vis.Vector3 = (1.0, 0.0, 0.0)
    'normal() -> tuple[float, float, float]\n\nThe normal vector of the slicing plane. Its interpretation depends on whether :py:attr:`miller` mode is enabled or not. \n\nDefault: ``(1.0, 0.0, 0.0)``'
    only_selected: bool = False
    'only_selected() -> bool\n\nControls whether the modifier should act only on currently selected data elements (e.g. selected particles).\n\nDefault: ``False``'
    operate_on: MutableSet[str] = {'particles', 'surfaces', 'voxels', 'lines', 'dislocations'}
    "operate_on() -> collections.abc.MutableSet[str]\n\nA set of strings specifying the kinds of data elements this modifier should operate on. By default the set contains all data element types supported by the modifier. \n\nDefault: ``{'particles', 'dislocations', 'lines', 'surfaces', 'vectors', 'voxels'}``"
    plane_vis: ovito.vis.TriangleMeshVis = ovito.vis.TriangleMeshVis(highlight_edges=True, transparency=0.5, title='Plane')
    'plane_vis() -> ovito.vis.TriangleMeshVis\n\nThe :py:class:`TriangleMeshVis` controlling the visual appearance of the cutting plane in rendered images. The visual element is only used when :py:attr:`render_plane` has been set to ``True`` to visualize the mathematical plane of the modifier as a visible polygon.'
    render_plane: bool = False
    'render_plane() -> bool\n\nControls whether the modifier should produce renderable geometry to visualize the cutting plane. The visual appearance of the plane can be adjusted through the :py:attr:`plane_vis` element. \n\nDefault: ``False``'
    select: bool = False
    'select() -> bool\n\nIf ``True``, the modifier selects data elements instead of deleting them.\n\nDefault: ``False``'
    slab_width: float = 0.0
    'slab_width() -> float\n\nThe thickness of the slab to cut (in simulation units of length). If zero, the modifier cuts away everything on one side of the cutting plane.\n\nDefault: ``0.0``'

@dataclass(kw_only=True)
class SmoothTrajectoryModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

This modifier smoothens the particle motion by averaging the particle positions from several successive snapshots of a simulation trajectory.
It can be used to create smooth-looking animations from a relatively coarse sequence of simulation snapshots
and reduce fluctuations or thermal vibrations in particle trajectories.
See also the corresponding user manual page for this modifier.

Note: Make sure you insert this modifier at the beginning of a data pipeline, in particular before any other modifiers that 
delete or filter particles, because needs to see the complete particle system in order to perform the trajectory smoothing."""
    minimum_image_convention: bool = True
    'minimum_image_convention() -> bool\n\nIf this option is set, the modifier will automatically detect when particles cross a simulation box boundary in between two successive simulation frames and computes the unwrapped displacements correctly. You should leave this option activated unless the input particle coordinates are already in unwrapped form. \n\nDefault: ``True``'
    window_size: int = 1
    'window_size() -> int\n\nControls how many input animation frames to take into account when calculating the time-averaged particle coordinates. The modifier uses a sliding time window of the given size that is centered around the current animation time. A window size of 1 disables smoothing. \n\nDefault: ``1``'

@dataclass(kw_only=True)
class SpatialBinningModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

The modifier outputs a 1-, 2- or 3-dimensional grid of field values. See also the corresponding
user manual page for more information.
This grid contains coarse-grained information of the :py:class:`Particles` or
a :py:class:`DislocationNetwork` found in the :py:class:`DataCollection`:

* Particles: This modifier maps the particles to the binning grid and applies a reduction operation
  (mean, sum, min, max, etc.) to some property of the particles located within each spatial bin.
  You have to specify the input particle property which should be transferred to the grid cell by
  setting the modifier's :py:attr:`property` field to the name of an existing particle property.

* DislocationNetwork: In this mode, the modifier calculates the local dislocation density and Nye tensor in each bin.
  The dislocation density is expressed in units of 1/length\\ :sup:`2`, depending on the unit system used in the
  input simulation.

The :py:attr:`direction` parameter of the modifier selects the dimensionality of the bin grid and the
spatial direction(s) along which the grid cells are oriented. The :py:attr:`bin_count` parameter controls the resolution
of the generated grid, which always covers the entire :py:class:`SimulationCell`.

If the dimensionality of the selected output grid is one-dimensional, then the modifier will generate a :py:class:`DataTable` object
containing the computed bin values along the selected spatial :py:attr:`direction`. The data table may later be retrieved from the
:py:attr:`tables` dictionary of the pipeline's output :py:class:`DataCollection`
under the lookup key ``'binning'``.

If the dimensionality of the selected output grid is two- or three-dimensional, then the modifier will generate a :py:class:`VoxelGrid` object.
The generated voxel grid may later be retrieved from the :py:attr:`grids` dictionary of the pipeline's output
:py:class:`DataCollection` under the lookup key ``'binning'``.

Examples

The following example shows how to calculate the gradient of the ``Velocity.X`` field along the spatial Z-axis.
For this, we use the :py:attr:`first_derivative` option, which lets the modifier compute the first derivative of the
data points using a finite differences scheme.

```python
  pipeline.modifiers.append(SpatialBinningModifier(
      property = 'Velocity.X',
      direction = SpatialBinningModifier.Direction.Z, 
      bin_count = 80,
      reduction_operation = SpatialBinningModifier.Operation.Mean,
      first_derivative = True
  ))
```

The output :py:class:`DataTable` can then be exported to an output file using OVITO's :py:func:`export_file`
function or printed to the console:

```python
  export_file(pipeline, 'output/velocity_profile.txt', 'txt/table', key='binning')
  data = pipeline.compute()
  print(data.tables['binning'].xy())
```

The following example demonstrates how to compute a three-dimensional :py:class:`VoxelGrid` of the particle number density.
Since the :py:class:`SpatialBinningModifier` always requires some input particle property, we first employ the
:py:class:`ComputePropertyModifier` to give all particles a new property with the uniform value 1, which can then
serve as input property for the binning:

```python
  pipeline.modifiers.append(ComputePropertyModifier(expressions=['1'], output_property='Unity'))
  
  pipeline.modifiers.append(SpatialBinningModifier(
      property = 'Unity',
      direction = SpatialBinningModifier.Direction.XYZ, 
      bin_count = (100, 100, 100),
      reduction_operation = SpatialBinningModifier.Operation.SumVol
  ))
```

The resulting :py:class:`VoxelGrid` may now be exported to an output file, for example in the VTK format:

```python
  export_file(pipeline, 'output/density.vtk', 'vtk/grid', key='binning')
```

Or it can serve as input for subsequent modifiers in your data pipeline, for example the :py:class:`CreateIsosurfaceModifier`."""

    class Direction(enum.Enum):
        """"""
        X = enum.auto()
        Y = enum.auto()
        Z = enum.auto()
        XY = enum.auto()
        XZ = enum.auto()
        YZ = enum.auto()
        XYZ = enum.auto()

    class Operation(enum.Enum):
        """"""
        Mean = enum.auto()
        Sum = enum.auto()
        SumVol = enum.auto()
        Min = enum.auto()
        Max = enum.auto()
    bin_count: Tuple[int, ...] = (50, 50, 20)
    "bin_count() -> Union[int, tuple[int, int], tuple[int, int, int]]\n\nSpecifies the number of bin cells to generate along each axis of the binning grid. You should assign a tuple containing one, two, or three positive integers to this parameter field, depending on the grid's dimensionality set by the :py:attr:`direction` parameter. \n\nNote that the entries in the tuple specify the number of bins along the grid's first, second, and third dimension and not along the spatial axes. For example, if the binning :py:attr:`direction` is ``Direction.YZ``, setting :py:attr:`bin_count` to ``(100, 50)`` will let the modifier generate a two-dimensional grid with 100 bins along the second simulation cell vector (spatial y-axis) and 50 bins along the third cell vector (z-axis) of the 3-dimensional simulation box. Examples:: \n\n   SpatialBinningModifier(direction=SpatialBinningModifier.Direction.Z, bin_count=100)\n   SpatialBinningModifier(direction=SpatialBinningModifier.Direction.XZ, bin_count=(40, 80))\n   SpatialBinningModifier(direction=SpatialBinningModifier.Direction.XYZ, bin_count=(80, 80, 80))\n\n\nIf you assign just a single number, or a tuple with fewer entries than required for the selected grid dimensionality, the bin count in any of the remaining grid dimensions will be implicitly set to 1. \n\nDefault: ``(50, 50, 20)``"
    direction: SpatialBinningModifier.Direction = Direction.X
    'direction() ->  SpatialBinningModifier.Direction\n\nSelects the alignment of the bins and the dimensionality of the grid. Possible values:\n\n   * ``SpatialBinningModifier.Direction.X``\n   * ``SpatialBinningModifier.Direction.Y``\n   * ``SpatialBinningModifier.Direction.Z``\n   * ``SpatialBinningModifier.Direction.XY``\n   * ``SpatialBinningModifier.Direction.XZ``\n   * ``SpatialBinningModifier.Direction.YZ``\n   * ``SpatialBinningModifier.Direction.XYZ``\n\nFor modes ``X``, ``Y``, and ``Z``, the modifier will generate a one-dimensional grid with bins aligned perpendicular to the selected simulation cell vector. For modes ``XY``, ``XZ``, and ``YZ``, the modifier will generate a two-dimensional grid with bins aligned perpendicular to both selected simulation cell vectors (i.e. parallel to the third vector). In the last case (``XYZ``), the modifier generates a three-dimensional voxel grid. \n\nDefault: ``SpatialBinningModifier.Direction.X``'
    first_derivative: bool = False
    'first_derivative() ->  bool\n\nSet this to ``True`` to let the modifier numerically compute the first derivative of the binned data points using a finite differences approximation. This works only if binning is performed in a single :py:attr:`direction` (``Direction.X``, ``Direction.Y`` or ``Direction.Z``). For two- or three-dimensional binning modes, this option is ignored. \n\nDefault: ``False``'
    grid_vis: ovito.vis.VoxelGridVis = ovito.vis.VoxelGridVis(title='Binning grid')
    'grid_vis() ->  ovito.vis.VoxelGridVis\n\nThe :py:class:`VoxelGridVis` element controlling the visual appearance of the 2d or 3d grid generated by the modifier in rendered images. The visual element is ignored if :py:attr:`direction` parameter is set to a 1-dimensional binning mode.'
    only_selected: bool = False
    'only_selected() ->  bool\n\nThis option lets the modifier take into account only the currently selected particles. Unselected particles will be excluded from the binning process. You can use this option to restrict the calculation to a subset of particles. \n\nDefault: ``False``'
    property: str = ''
    'property() ->  str\n\nThe name of the input particle :py:class:`Property` which the reduction operation should be applied to. This can be the name of one of the standard particle properties or of a user-defined particle property. \n\nIf no input input property is specified, then the modifier uses unity as uniform input value for all particles. Use this to compute the number density of particles in each bin, for example, with :py:attr:`reduction_operation` set to ``SumVol``. \n\nIf the input is a vector property, a component name may be appended after a dot, e.g. ``"Velocity.X"``, to perform the reduction operation only on that specific vector component. The output will then be a scalar field. Otherwise, the reduction operation is applied to all vector components independently and the output will be a vector field or vector-valued function. \n\nDefault: ``\'\'``'
    operate_on: str = 'particles'
    "operate_on() -> str\n\nSelects the data element this modifier will operate on. Valid values are: \n\n  * ``'particles'``\n  * ``'dislocations'``\n\n\nDefault: ``'particles'``"
    reduction_operation: SpatialBinningModifier.Operation = Operation.Mean
    'reduction_operation() ->  SpatialBinningModifier.Operation\n\nSelects the reduction operation to be carried out. Supported parameter values are:\n\n   * ``SpatialBinningModifier.Operation.Mean``\n   * ``SpatialBinningModifier.Operation.Sum``\n   * ``SpatialBinningModifier.Operation.SumVol``\n   * ``SpatialBinningModifier.Operation.Min``\n   * ``SpatialBinningModifier.Operation.Max``\n\nThe operation ``SumVol`` first computes the sum and then divides the result by the volume of the respective bin. It is intended to compute the number density of particles or pressure/stress within each bin from the per-atom virial. \n\nDefault: ``SpatialBinningModifier.Operation.Mean``'

@dataclass(kw_only=True)
class SpatialCorrelationFunctionModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

This modifier calculates the spatial correlation function between two particle properties. See also the corresponding user manual page for this modifier. 

The algorithm uses the FFT to compute the convolution. It then computes a radial average in reciprocal and real space. This gives the correlation function up to half of the cell size. The modifier can additionally compute the short-ranged part of the correlation function from a direct summation over neighbors.

Usage example:

```python
  from ovito.modifiers import SpatialCorrelationFunctionModifier
  from ovito.io import import_file, export_file
  
  pipeline = import_file("input/simulation.dump")
  
  mod = SpatialCorrelationFunctionModifier(property1='Position.X', property2='peatom')
  pipeline.modifiers.append(mod)
  data = pipeline.compute()
  
  # Export RDF and correlation functions to text files:
  C = data.tables['correlation-real-space']
  rdf = data.tables['correlation-real-space-rdf']
  export_file(C, 'output/real_correlation_function.txt', 'txt/table')
  export_file(rdf, 'output/rdf.txt', 'txt/table')
  
  # Compute normalized correlation function (by co-variance):
  C_norm = C.xy()
  mean1 = data.attributes['CorrelationFunction.mean1']
  mean2 = data.attributes['CorrelationFunction.mean2']
  covariance = data.attributes['CorrelationFunction.covariance']
  C_norm[:,1] = (C_norm[:,1] - mean1*mean2) / (covariance - mean1*mean2)
  import numpy
  numpy.savetxt('output/normalized_real_correlation_function.txt', C_norm)
```"""
    apply_window: bool = True
    'apply_window() -> bool\n\nThis flag controls whether non-periodic directions have a Hann window applied to them. Applying a window function is necessary to remove spurious oscillations and power-law scaling of the (implicit) rectangular window of the non-periodic domain. \n\nDefault: ``True``'
    direct_summation: bool = False
    'direct_summation() -> bool\n\nFlag controlling whether the real-space correlation plot will show the result of a direct calculation of the correlation function, obtained by summing over neighbors. \n\nDefault: ``False``'
    grid_spacing: float = 3.0
    'grid_spacing() -> float\n\nControls the approximate size of the FFT grid cell. The actual size is determined by the distance of the simulation cell faces which must contain an integer number of grid cells. \n\nDefault: ``3.0``'
    neighbor_bins: int = 50
    'neighbor_bins() -> int\n\nThis integer value controls the number of bins for the direct calculation of the real-space correlation function. \n\nDefault: ``50``'
    neighbor_cutoff: float = 5.0
    'neighbor_cutoff() -> float\n\nThis parameter determines the cutoff of the direct calculation of the real-space correlation function. \n\nDefault: ``5.0``'
    property1: str = ''
    'property1() -> str\n\nThe name of the first input particle property for which to compute the correlation, P1. For vector properties a component name must be appended in the string, e.g. ``"Velocity.X"``. \n\nDefault: ``\'\'``'
    property2: str = ''
    "property2() -> str\n\nThe name of the second particle property for which to compute the correlation, P2. If this is the same as :py:attr:`property1`, then the modifier will compute the autocorrelation. \n\nDefault: ``''``"

@dataclass(kw_only=True)
class TimeAveragingModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

This modifier computes the time average of one or more time-dependent input quantities by sampling them over all frames of the loaded trajectory.
The following types of input quantities are supported:

  * time-dependent global attributes
  * :py:class:`DataTable` objects, e.g. histograms or distribution functions computed by other analysis modifiers in the pipeline
  * element-wise property values, e.g. particle coordinates or other particle properties
  * :py:class:`SimulationCell` (the cell shape)

The selected input quantities to be averaged must be produced in the upstream data pipeline or loaded from the input trajectory file.
See also the corresponding user manual page for more information on this modifier.

Averaging scalar quantities:

The following code example demonstrates how to use the :py:class:`TimeAveragingModifier` to calculate the mean value of
a scalar quantity that changes with simulation time. In this example, the :py:class:`CreateBondsModifier` is
added to the pipeline to dynamically create bonds between the moving particles (thus, the number of created bonds will change with time).
The modifier reports the number of bonds that have been created in the current simulation frame as an output attribute named ``'CreateBonds.num_bonds'``.
This varying attribute is then selected as the input quantity to be averaged by the following :py:class:`TimeAveragingModifier`,
which in turn outputs the time-averaged value it computes as a new attribute named ``'CreateBonds.num_bonds (Average)'``.
This value can be retrieved from the :py:class:`DataCollection` produced by `Pipeline.compute()`.
Note that is does not matter at which simulation time we evaluate the data pipeline -- the mean number of bonds is a static quantity and remains
constant over the entire trajectory.

```python
  pipeline = import_file('input/trajectory.dump')
  pipeline.modifiers.append(CreateBondsModifier(cutoff=2.9))
  pipeline.modifiers.append(TimeAveragingModifier(operate_on='attribute:CreateBonds.num_bonds'))
  data = pipeline.compute()
  print(data.attributes['CreateBonds.num_bonds (Average)'])
```

Averaging data tables:

Some analysis modifiers of OVITO output their results as :py:class:`DataTable` objects, which can represent
histograms or distributions functions. For example, the :py:class:`CoordinationAnalysisModifier` computes
the radial pair distribution function (RDF) of a particle system for the current simulation time and outputs it as a data table
named ``'coordination-rdf'``. We can have the :py:class:`TimeAveragingModifier` compute the time average of this varying
distribution function over all frames of the loaded trajectory:

```python
  pipeline = import_file('input/trajectory.dump')
  pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff=6.0, number_of_bins=50))
  pipeline.modifiers.append(TimeAveragingModifier(operate_on='table:coordination-rdf'))
  data = pipeline.compute()
  print(data.tables['coordination-rdf[average]'].xy())
```

Note that, in case of data tables, the output table's :py:attr:`identifier` is given the suffix ``[average]``
by the :py:class:`TimeAveragingModifier`.

An average of a :py:class:`DataTable` can only be computed if the table's x-axis :py:attr:`interval` does not vary with time and
the :py:attr:`x`-coordinates of the data points remain fixed. That's because the modifier simply averages the time-varying y-coordinate of each data point.
Thus, to generate a histogram :py:class:`DataTable` using the :py:class:`HistogramModifier` that is suitable for time averaging,
you will have to activate the :py:attr:`HistogramModifier.fix_xrange` option.

Averaging properties:

The modifier can time-average properties that belong to a :py:class:`Particles` object, :py:class:`VoxelGrid`,
or other type of :py:class:`PropertyContainer`. When setting the :py:attr:`operate_on` field, you need to specify
the container's identifier and the name of the property to average:

```python
  pipeline.modifiers.append(TimeAveragingModifier(operate_on='property:particles/Coordination'))
  data = pipeline.compute()
  print(data.particles['Coordination Average'][...])
```

Note that the output :py:attr:`Property` is given the suffix ``'Average'`` by the :py:class:`TimeAveragingModifier`."""
    interval: Optional[Tuple[int, int]] = None
    'interval() -> Optional[tuple[int, int]]\n\nThe animation frame interval over which the input trajectory is sampled to compute the average. You can set this to a pair of integers specifying the first and the last frame of the averaging interval; or assign ``None`` to let the modifier compute the average over the entire trajectory. \n\nFor example, to restrict the average to the second half of the loaded simulation trajectory, you can specify the frame interval based on the :py:attr:`num_frames` value: \n\n```python\n  modifier.interval = (pipeline.num_frames//2, pipeline.num_frames-1)\n```\n\nDefault: ``None``'
    operate_on: Union[str, Sequence[str]] = ''
    "operate_on() -> Union[str, collections.abc.Sequence[str]]\n\nSelects the input quantity to be averaged by this modifier. Supported values for this field are: \n\n  * ``'attribute:<NAME>'``\n  * ``'table:<ID>'``\n  * ``'property:<CONTAINER>/<PROPERTY>'``\n  * ``'cell'``\n\n\nHere, ``<NAME>`` refers to the name of a global attribute to be averaged. ``<ID>`` is the :py:attr:`identifier` of the :py:class:`DataTable` to be averaged, ``<CONTAINER>`` is the :py:attr:`identifier` of a :py:class:`PropertyContainer` and ``<PROPERTY>`` the :py:attr:`name` of the :py:class:`Property` in that container, which should be averaged. Note that the :py:class:`Particles` property container has the standard identifier ``'particles'``. Its :py:class:`Bonds` child container is referenced by the hierarchical identifier ``'particles/bonds'``. \n\nFurthermore, it is possible to let the modifier average multiple input quantities in one pass by specifying a list of input references of the kind described above. The following code demonstrates how to compute the time averages of two global attributes: \n\n```python\n  pipeline.modifiers.append(TimeAveragingModifier(\n      operate_on = ('attribute:CommonNeighborAnalysis.counts.FCC',\n                    'attribute:CommonNeighborAnalysis.counts.BCC')\n  ))\n  data = pipeline.compute()\n  print(data.attributes['CommonNeighborAnalysis.counts.FCC (Average)'])\n  print(data.attributes['CommonNeighborAnalysis.counts.BCC (Average)'])\n```\n\nDefault: ``''``"
    overwrite: bool = False
    'overwrite() -> bool\n\nIf set to false, the averaged values computed by the modifier are output as a new quantities with the suffix "Average" appended to keep the original and averaged values separate. If set to true, the modifier will overwrite the current time-dependent values in the trajectory with the static averaged values. \n\nDefault: ``False``'
    sampling_frequency: int = 1
    'sampling_frequency() -> int\n\nAnimation step interval at which the quantity to be averaged is sampled from the input trajectory. You can set this to a larger value to perform a coarser sampling and reduce the total number of trajectory frames that need to be loaded into memory. \n\nDefault: ``1``'

@dataclass(kw_only=True)
class TimeSeriesModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

This modifier samples a varying input quantity computed by the data pipeline over all frames of the loaded trajectory
to produce a time series, which can then be plotted to visualize the time evolution of the quantity.
See also the corresponding user manual page for more information on this modifier. 

The modifier takes one or more global attributes as input, samples their values over the simulation trajectory, 
and outputs the generated time series as a new :py:class:`DataTable` with the identifier ``time-series``. 

```python
  pipeline = import_file('input/trajectory.dump')
  pipeline.modifiers.append(ConstructSurfaceModifier(radius = 1.9))
  pipeline.modifiers.append(TimeSeriesModifier(operate_on = 'ConstructSurfaceMesh.surface_area'))
  data = pipeline.compute()
  print(data.tables['time-series'].xy())
```

This code samples the value of the global attribute ``ConstructSurfaceMesh.surface_area``, which is computed by the :py:class:`ConstructSurfaceModifier`
for each frame of the loaded simulation trajectory. The result is a static :py:class:`DataTable` with the number of rows equal to the number
of trajectory frames, which is inserted by the :py:class:`TimeSeriesModifier` into the pipeline's output data collection.

Note that you do not need the :py:class:`TimeSeriesModifier` in order to output a global attribute to a text file as a function of time.
This can be accomplished simply by calling the :py:func:`export_file` function, which can automatically sample an attribute's value 
and write it to the output file one line per frame:

```python
  pipeline = import_file('input/trajectory.dump')
  pipeline.modifiers.append(ConstructSurfaceModifier(radius = 1.9))
  export_file(pipeline, 'output/surface_area.txt',
      format='txt/attr',
      columns=["SourceFrame", "ConstructSurfaceMesh.surface_area"],
      multiple_frames=True)
```"""
    interval: Optional[Tuple[int, int]] = None
    'interval() -> Optional[tuple[int, int]]\n\nThe interval of animation frames over which the input trajectory is sampled to generate the time series. You can set this to a pair of integers specifying the first and the last frame of the sampling interval; or assign ``None`` to let the modifier generate the time series over the entire trajectory. \n\nFor example, to restrict the time series to the first half of the loaded simulation trajectory, you can specify the frame interval based on the :py:attr:`num_frames` value: \n\n```python\n  modifier.interval = (0, pipeline.num_frames//2)\n```\n\nDefault: ``None``'
    operate_on: Union[str, Sequence[str]] = ''
    'operate_on() -> Union[str, collections.abc.Sequence[str]]\n\nSpecifies the name of the input global attribute to be sampled by the modifier. The attribute must be generated by the trajectory file reader or dynamically computed by a modifier in the upstream data pipeline. \n\nYou can also specify several global attributes as a Python tuple of strings. Then the modifier will simultaneously sample each of the input attributes and produce a :py:class:`DataTable` with a vector data column containing the set of time series. \n\n```python\n  pipeline.modifiers.append(TimeSeriesModifier(\n      operate_on = (\'ConstructSurfaceMesh.surface_area\',\n                    \'ConstructSurfaceMesh.filled_volume\')))\n  data = pipeline.compute()\n  series = data.tables[\'time-series\'].y\n  print("Surface area:", series[:,0])\n  print("Solid volume:", series[:,1])\n```\n\nDefault: ``\'\'``'
    sampling_frequency: int = 1
    'sampling_frequency() -> int\n\nAnimation interval at which the attribute(s) is sampled from the input trajectory. You can set this to a larger value to perform a coarser sampling and reduce the total number of trajectory frames that need to be loaded/computed. \n\nDefault: ``1``'
    time_attribute: str = ''
    "time_attribute() -> str\n\nThe name of a global attribute that should be queried by the modifier to obtain the time-axis coordinates of the data samples. If set to an empty string (the default), the modifier uses the animation frame number as time axis. You can alternatively tell the modifier to use the ``Timestep`` or ``Time`` global attributes, which are created by some file readers based on the information found in the input trajectory, in order to plot the selected input attribute as a function of  MD timestep number or physical simulation time. \n\nDefault: ``''``"

@dataclass(kw_only=True)
class UnwrapTrajectoriesModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

This modifier determines when particles cross through the periodic boundaries of the simulation cell and unwraps the particle coordinates in order to make the trajectories continuous. As a result of this operation, particle trajectories will no longer fold back into the simulation cell and instead lead outside the cell. 

Note that, to unwrap the particle coordinates, the modifier may have to step through all frames of the input simulation trajectory to detect jumps in the trajectories. This will not be necessary, however, if the ``Periodic Image`` particle property has been loaded from the input simulation file, because then the modifier can directly use this information to unwrap the particle coordinates. 

```python
  from ovito.io import import_file
  from ovito.modifiers import UnwrapTrajectoriesModifier
  
  # Load a simulation trajectory:
  pipeline = import_file('input/simulation.*.dump')
  
  # Insert the unwrap modifier into the pipeline.
  # Note that the modifier should typically precede any other modifiers in the pipeline.
  pipeline.modifiers.append(UnwrapTrajectoriesModifier())
  
  # For demonstration purposes, request last frame of the trajectory.
  # The returned DataCollection will contain modified particle coordinates, which are computed
  # by the unwrap modifier by tracing the trajectories of the particles and unfolding
  # them whenever crossings of the periodic cell boundaries are detected.
  data = pipeline.compute(pipeline.num_frames - 1)
```"""
    pass

@dataclass(kw_only=True)
class VoroTopModifier(StructureIdentificationModifier):
    """Base: :py:class:`ovito.modifiers.StructureIdentificationModifier`

"This modifier uses the Voronoi cell topology of particles to characterize their local structural environments
"[`Lazar, Han, Srolovitz, PNAS 112:43 (2015) <http://dx.doi.org/10.1073/pnas.1505788112>`__].
See the corresponding user manual page
for more information on this modifier. Note that this modifier inherits several important parameter fields 
from its :py:class:`StructureIdentificationModifier` base class.

The Voronoi cell of a particle is the region of space closer to it than to any other particle. 
The topology of the Voronoi cell is the manner in which its faces are connected, and describes 
the manner in which a particle's neighbors are arranged.  The topology of a Voronoi cell can be 
completely described in a vector of integers called a *Weinberg vector* 
[`Weinberg, IEEE Trans. Circuit Theory 13:2 (1966) <http://dx.doi.org/10.1109/TCT.1966.1082573>`__]. 

This modifier requires loading a *filter*, which specifies structure types and associated 
Weinberg vectors.  Filters for several common structures can be obtained from the 
`VoroTop <https://www.vorotop.org/download.html>`__ website. 
The modifier calculates the Voronoi cell topology of each particle, uses the provided 
filter to determine the structure type, and stores the results in the ``Structure Type`` particle property. 
This allows the user to subsequently select particles  of a certain structural type, e.g. by using the 
:py:class:`SelectTypeModifier`. 

This method is well-suited for analyzing finite-temperature systems, including those heated to 
their bulk melting temperatures. This robust behavior relieves the need to quench a sample 
(such as by energy minimization) prior to analysis. 
Further information about the Voronoi topology approach for local structure analysis, as well 
as additional filters, can be found on the `VoroTop webpage <https://www.vorotop.org/>`__."""
    filter_file: str = ''
    "filter_file() -> str\n\nPath to the filter definition file used by the modifier. Filters files are available from the `VoroTop <https://www.vorotop.org/download.html>`__ website. \n\nDefault: ``''``"
    use_radii: bool = False
    'use_radii() -> bool\n\nIf ``True``, the modifier computes the poly-disperse Voronoi tessellation, which takes into account the radii of particles. Otherwise a mono-disperse Voronoi tessellation is computed, which is independent of the particle sizes. \n\nDefault: ``False``'

@dataclass(kw_only=True)
class VoronoiAnalysisModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

Computes the atomic volumes and coordination numbers using a Voronoi tessellation of the particle system.
See the corresponding user manual page for more information.
See :ref:`example_compute_voronoi_indices` for a code example demonstrating the use of this modifier.

Inputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Position``
      - The coordinates of the input particles.
    * - ``Radius``
      - Per-particle radii are used if :py:attr:`use_radii` is set to ``True``.
    * - ``Particle Type``
      - Per-type radii are used if :py:attr:`use_radii` is set to ``True`` and ``Radius`` property is not present.
    * - ``Selection``
      - The selection state of the input particles. Only needed if :py:attr:`only_selected` is set to ``True``.

Outputs:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Particle properties
      -
    * - ``Atomic Volume``
      - The computed volume of each particle's Voronoi polyhedron.
    * - ``Coordination``
      - The number of faces of each particle's Voronoi polyhedron.
    * - ``Voronoi Index``
      - The index vector of each Voronoi polyhedron. Only computed if :py:attr:`compute_indices` is set to ``True``.
    * - ``Max Face Order``
      - The maximum number of edges in any face of a particle's Voronoi polyhedron. Only if :py:attr:`compute_indices` is set to ``True``.
    * - ``Cavity Radius``
      - Distance from a center particle to the farthest vertex of its Voronoi polyhedron.

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Global attributes
      -
    * - ``Voronoi.max_face_order``
      - Indicates the maximum number of edges of any face in the computed Voronoi tessellation (ignoring edges and faces that fall below the area/length thresholds).

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Bond properties
      -
    * - ``Topology``
      - The connectivity information of newly created :py:class:`Bonds` (one bond for each Voronoi face). Only if :py:attr:`generate_bonds` is set to ``True``."""
    bonds_vis: ovito.vis.BondsVis = ovito.vis.BondsVis()
    'bonds_vis() -> ovito.vis.BondsVis\n\nThe :py:class:`BondsVis` object controlling the visual appearance of the bonds generated by the modifier if :py:attr:`generate_bonds` is set to ``True``.'
    compute_indices: bool = False
    'compute_indices() -> bool\n\nIf ``True``, the modifier calculates the Voronoi indices of particles. The modifier stores the computed indices in a vector particle property named ``Voronoi Index``. The *i*-th component of this property will contain the number of faces of the Voronoi cell that have *i* edges. Thus, the first two components of the per-particle vector will always be zero, because the minimum number of edges a polygon can have is three. \n\nDefault: ``False``'
    edge_threshold: float = 0.0
    'edge_threshold() -> float\n\nSpecifies the minimum length an edge must have to be considered in the Voronoi index calculation. Edges that are shorter than this threshold will be ignored when counting the number of edges of a Voronoi face. The threshold parameter is an absolute value in units of length of your input data. \n\nDefault: ``0.0``'
    face_threshold: float = 0.0
    'face_threshold() -> float\n\nSpecifies a minimum area for individual Voronoi faces in terms of an absolute area. The algorithm will ignore any face of a Voronoi polyhedron with an area smaller than this threshold when computing the coordination number and the Voronoi index of a particle. The threshold parameter is an absolute area given in units of length squared (in whatever units your input data is given). \n\nNote that this absolute area threshold and the :py:attr:`relative_face_threshold` are applied simultaneously. \n\nDefault: ``0.0``'
    generate_bonds: bool = False
    "generate_bonds() -> bool\n\nControls whether the modifier outputs the nearest neighbor bonds. The modifier will generate a bond for every pair of adjacent atoms that share a face of the Voronoi tessellation. No bond will be created if the face's area is below the :py:attr:`face_threshold` or if the face has less than three edges that are longer than the :py:attr:`edge_threshold`.\n\nDefault: ``False``"
    generate_polyhedra: bool = False
    'generate_polyhedra() -> bool\n\nControls whether the modifier outputs the computed Voronoi cells as a polyhedral :py:class:`SurfaceMesh` object. \n\nDefault: ``False``'
    mesh_vis: ovito.vis.SurfaceMeshVis = ovito.vis.SurfaceMeshVis(show_cap=False, smooth_shading=False, surface_transparency=0.25, highlight_edges=True, title='Voronoi polyhedra')
    'mesh_vis() -> ovito.vis.SurfaceMeshVis\n\nThe :py:class:`SurfaceMeshVis` object controlling the visual appearance of the polyhedral mesh generated by the modifier if :py:attr:`generate_polyhedra` is set to ``True``.'
    only_selected: bool = False
    'only_selected() -> bool\n\nLets the modifier perform the analysis only for selected particles. Particles that are currently not selected will be treated as if they did not exist.\n\nDefault: ``False``'
    relative_face_threshold: float = 0.0
    'relative_face_threshold() -> float\n\nSpecifies a minimum area for Voronoi faces in terms of a fraction of total area of the Voronoi polyhedron surface. The algorithm will ignore any face of a Voronoi polyhedron with an area smaller than this threshold when computing the coordination number and the Voronoi index of particles. The threshold parameter is specified as a fraction of the total surface area of the Voronoi polyhedron the faces belong to. For example, a threshold value of 0.01 would remove those faces from the analysis with an area less than 1% of the total area of the polyhedron surface. \n\nNote that this relative threshold and the absolute :py:attr:`face_threshold` are applied simultaneously. \n\nDefault: ``0.0``'
    use_radii: bool = False
    'use_radii() -> bool\n\nIf ``True``, the modifier computes the poly-disperse Voronoi tessellation, which takes into account the radii of particles. Otherwise a mono-disperse Voronoi tessellation is computed, which is independent of the particle sizes. \n\nDefault: ``False``'

@dataclass(kw_only=True)
class WignerSeitzAnalysisModifier(ovito.pipeline.ReferenceConfigurationModifier):
    """Base: :py:class:`ovito.pipeline.ReferenceConfigurationModifier`

Performs the Wigner-Seitz cell analysis to identify point defects in a crystal.
See the corresponding user manual page for more information.

Defects are identified with respect to a perfect reference crystal configuration.
By default, frame 0 of the current simulation sequence is used as reference configuration.
The modifier inherits from the :py:class:`ReferenceConfigurationModifier` class, which provides
further settings that control the definition of the reference configuration.

Outputs:

.. list-table::
    :widths: 20 80
    :header-rows: 1

    * - Particle properties
      -
    * - ``Occupancy``
      - The computed site occupation numbers, one for each site in the reference configuration, if :py:attr:`output_displaced` is set to false.
        Otherwise, total number of atoms occupying the same reference site as the atom from the displaced configuration.
    * - ``Site Index``
      - Zero-based index of the reference site to which the atom from the displaced configuration has been assigned. Only available if :py:attr:`output_displaced` is set to true.
    * - ``Site Identifier``
      - Unique identifier of the reference site to which the atom from the displaced configuration has been assigned. Only available if :py:attr:`output_displaced` is set to true and if reference sites have an `Particle Identifier` property.
    * - ``Site Type``
      - Type of the reference site to which the atom from the displaced configuration has been assigned. Only available if :py:attr:`output_displaced` is set to true.

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Global attributes
      -
    * - ``WignerSeitz.vacancy_count``
      - The total number of vacant sites (having ``Occupancy`` == 0).
    * - ``WignerSeitz.interstitial_count``
      - The total number of of interstitial atoms. This is equal to the sum of occupancy numbers of all non-empty sites minus the number of non-empty sites.

Example:

The ``Occupancy`` particle property generated by the Wigner-Seitz algorithm allows you to select specific types of point defects, e.g.
antisites, using OVITO's selection tools. One option is to use the :py:class:`ExpressionSelectionModifier` to pick
sites having a certain occupancy. The following script exemplarily shows how a custom :py:class:`PythonModifier` function can be used to
select and count A-sites occupied by B-atoms in a binary system with two atom types (A=1 and B=2).

```python
  from ovito.io import *
  from ovito.data import *
  from ovito.modifiers import *
  from ovito.pipeline import *
  import numpy as np
  
  pipeline = import_file("input/simulation.*.dump")
  
  # Perform Wigner-Seitz analysis:
  ws = WignerSeitzAnalysisModifier(
      per_type_occupancies = True,
      affine_mapping = ReferenceConfigurationModifier.AffineMapping.ToReference)
  pipeline.modifiers.append(ws)
  
  # Define a modifier function that selects sites of type A=1 which
  # are occupied by exactly one atom of type B=2.
  def modify(frame, data):
  
      # Retrieve the two-dimensional array with the site occupancy numbers.
      # Use [...] to cast it to a Numpy array.
      occupancies = data.particles['Occupancy']
  
      # Get the site types as additional input:
      site_type = data.particles['Particle Type']
  
      # Calculate total occupancy of every site:
      total_occupancy = np.sum(occupancies, axis=1)
  
      # Set up a particle selection by creating the Selection property:
      selection = data.particles_.create_property('Selection')
  
      # Select A-sites occupied by exactly one B-atom (the second entry of the Occupancy
      # array must be 1, and all others 0). Note that the Occupancy array uses 0-based
      # indexing, while atom type IDs are typically 1-based.
      selection[...] = (site_type == 1) & (occupancies[:,1] == 1) & (total_occupancy == 1)
  
      # Additionally output the total number of antisites as a global attribute:
      data.attributes['Antisite_count'] = np.count_nonzero(selection)
  
  # Insert Python modifier into the data pipeline.
  pipeline.modifiers.append(modify)
  
  # Let OVITO do the computation and export the number of identified
  # antisites as a function of simulation time to a text file:
  export_file(pipeline, "output/antisites.txt", "txt/attr",
      columns = ['Timestep', 'Antisite_count'],
      multiple_frames = True)
  
  # Export the XYZ coordinates of just the antisites by removing all other atoms.
  pipeline.modifiers.append(InvertSelectionModifier())
  pipeline.modifiers.append(DeleteSelectedModifier())
  export_file(pipeline, "output/antisites.xyz", "xyz",
      columns = ['Position.X', 'Position.Y', 'Position.Z'],
      multiple_frames = True)
```"""
    output_displaced: bool = False
    'output_displaced() -> bool\n\nSpecifies whether the modifier should output the atoms of the current configuration or replace them with the sites from the reference configuration. \n\nBy default, the modifier throws away all atoms of the current configuration and outputs the atomic sites from the reference configuration instead. Thus, in this default mode, you will obtain information about how many atoms occupy each site from the reference configuration. If, however, you are more interested in visualizing the physical atoms that are currently occupying the sites (instead of the sites being occupied), then you should activate this modifier option. If set to true, the modifier will maintain the input atoms from the current configuration. The ``Occupancy`` property generated by the modifier will now pertain to the atoms instead of the sites, with the following new meaning: The occupancy number now counts how many atoms in total are occupying the same site as the atom the property refers to does. Furthermore, the modifier will in this mode output another property named ``Site Type``, which reports for each atom the type of the reference site it was assigned to by the W-S algorithm. \n\nDefault: ``False``'
    per_type_occupancies: bool = False
    'per_type_occupancies() -> bool\n\nA flag controlling whether the modifier should compute occupancy numbers on a per-particle-type basis. \n\nIf false, only the total occupancy number is computed for each reference site, which counts the number of particles that occupy the site irrespective of their types. If true, then the ``Occupancy`` property computed by the modifier becomes a vector property with *N* components, where *N* is the number of particle types defined in the system. Each component of the ``Occupancy`` property counts the number of particles of the corresponding type that occupy the site. For example, the property component ``Occupancy.1`` contains the number of particles of type 1 that occupy a site. \n\nDefault: ``False``'

@dataclass(kw_only=True)
class WrapPeriodicImagesModifier(ovito.pipeline.Modifier):
    """Base: :py:class:`ovito.pipeline.Modifier`

This modifier maps particles located outside of the simulation cell back into the cell by "wrapping" their coordinates 
around at the periodic boundaries of the :py:class:`SimulationCell`. 
See also the corresponding user manual page for this modifier. 
This modifier has no configurable parameters."""
    pass