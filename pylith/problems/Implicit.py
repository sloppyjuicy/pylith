#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
#

## @file pylith/problems/Implicit.py
##
## @brief Python Implicit object for solving equations using an
## implicit formulation.
##
## Factory: pde_formulation

from Formulation import Formulation

# Implicit class
class Implicit(Formulation):
  """
  Python Implicit object for solving equations using an implicit
  formulation.

  The formulation has the general form, [A(t+dt)] {u(t+dt)} = {b(t+dt)},
  where we want to solve for {u(t+dt)}. [A(t+dt)] generally
  depends on {u(t+dt)} as well as the current stresses and additional
  state variables.  For linear elastic or viscoelastic problems with
  constant time step size, A is a constant (after the elastic solution).
  {b(t+dt)} generally depends on the loads applied for time step t+dt
  (including the contributions to the internal force vector from
  displacement/velocity BC) as well as the internal force vector computed
  from the current stresses.

  Jacobian: A(t+dt)
  solution: u(t+dt)
  residual: bextern(t+dt) - bintern(t+dt)

  Factory: pde_formulation.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Formulation.Inventory):
    """
    Python object for managing Implicit facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Implicit facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="implicit"):
    """
    Constructor.
    """
    Formulation.__init__(self, name)
    self.outputField = {'name': "dispTBctpdt",
                        'label': "displacements"}
    return


  def elasticityIntegrator(self):
    """
    Get integrator for elastic material.
    """
    from pylith.feassemble.ElasticityImplicit import ElasticityImplicit
    return ElasticityImplicit()


  def initialize(self, mesh, materials, boundaryConditions,
                 interfaceConditions, dimension, dt):
    """
    Initialize problem for implicit time integration.
    """
    Formulation.initialize(self, mesh, materials, boundaryConditions,
                           interfaceConditions, dimension, dt)

    self._info.log("Creating other fields and matrices.")
    self.fields.addReal("dispTBctpdt")
    self.fields.addReal("residual")
    self.fields.copyLayout("solution")
    self.jacobian = mesh.createMatrix(self.fields.getReal("solution"))
    self.solver.initialize(mesh, self.fields.getReal("solution"))

    # Initial time step solves for total displacement field, not increment
    for constraint in self.constraints:
      constraint.useSolnIncr(False)
    for integrator in self.integrators:
      integrator.useSolnIncr(False)
    return


  def stableTimeStep(self):
    """
    Get stable time step for advancing forward in time.
    """
    self._info.log("WARNING: Implicit::stableTimeStep() not implemented.")
    from pyre.units.time import second
    dt = 0.0*second
    return dt
  

  def prestep(self, t, dt):
    """
    Hook for doing stuff before advancing time step.
    """
    self._info.log("Preparing for advanving time step.")

    # Set dispTBctpdt to the BC t time t+dt. Unconstrained DOF are
    # unaffected and will be equal to their values at time t.
    dispTBctpdt = self.fields.getReal("dispTBctpdt")
    for constraint in self.constraints:
      constraint.setField(t+dt, dispTBctpdt)

    needNewJacobian = False
    for integrator in self.integrators:
      if integrator.needNewJacobian():
        needNewJacobian = True
    if needNewJacobian:
      self._info.log("Reforming Jacobian of operator.")
      import pylith.utils.petsc as petsc
      petsc.mat_setzero(self.jacobian)
      for integrator in self.integrators:
        integrator.timeStep(dt)
        integrator.integrateJacobian(self.jacobian, t, self.fields)
      petsc.mat_assemble(self.jacobian)
    return


  def step(self, t, dt):
    """
    Advance to next time step.
    """
    self._info.log("Integrating residual term in operator.")
    residual = self.fields.getReal("residual")
    dispIncr = self.fields.getReal("solution")
    import pylith.topology.topology as bindings
    bindings.zeroRealSection(residual)
    bindings.zeroRealSection(dispIncr)
    for integrator in self.integrators:
      integrator.timeStep(dt)
      integrator.integrateResidual(residual, t, self.fields)

    #bindings.sectionView(residual, "residual")
    self._info.log("Solving equations.")
    self.solver.solve(dispIncr, self.jacobian, residual)
    #bindings.sectionView(self.fields.getReal("solution"), "solution")
    return


  def poststep(self, t):
    """
    Hook for doing stuff after advancing time step.
    """
    # After solving, dispTBctpdt contains the displacements at time t
    # for unconstrained DOF and displacements at time t+dt at
    # constrained DOF. We add in the displacement increments (only
    # nonzero at unconstrained DOF) so that after poststep(),
    # dispTBctpdt contains the displacement field at time t+dt.
    import pylith.topology.topology as bindings
    solution = self.fields.getReal("solution")
    disp = self.fields.getReal("dispTBctpdt")
    bindings.addRealSections(disp, disp, solution)
    #bindings.sectionView(solution, "solution")

    self._info.log("Updating integrators states.")
    for integrator in self.integrators:
      integrator.updateState(t, disp)

    # If finishing first time step, then switch from solving for total
    # displacements to solving for incremental displacements
    if 1 == self._istep:
      for constraint in self.constraints:
        constraint.useSolnIncr(True)
      for integrator in self.integrators:
        integrator.useSolnIncr(True)

    Formulation.poststep(self, t)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Formulation._configure(self)
    return


  def _solveElastic(self, mesh, materials, t, dt):
    """
    Solve for elastic solution.
    """
    self._info.log("Computing elastic solution.")

    self._info.log("Setting constraints.")
    solnField = self.fields.getReal("dispTBctpdt")
    import pylith.topology.topology as bindings
    bindings.zeroRealSection(solnField)
    for constraint in self.constraints:
      constraint.setField(t, solnField)

    self._info.log("Integrating Jacobian and residual of operator.")
    import pylith.utils.petsc as petsc
    petsc.mat_setzero(self.jacobian)
    residual = self.fields.getReal("residual")
    petsc.mat_setzero(self.jacobian)
    bindings.zeroRealSection(residual)
    for integrator in self.integrators:
      integrator.timeStep(dt)
      integrator.integrateJacobian(self.jacobian, t, self.fields)
      integrator.integrateResidual(residual, t, self.fields)
    import pylith.utils.petsc as petsc
    petsc.mat_assemble(self.jacobian)

    import pylith.topology.topology as bindings
    #petsc.mat_view(self.jacobian)
    #bindings.sectionView(residual, "residual")

    self._info.log("Solving equations.")
    self.solver.solve(solution, self.jacobian, residual)

    bindings.addRealSections(disp, disp, solution)
    #bindings.sectionView(solution, "solution")
    #bindings.sectionView(self.fields.getReal("dispIncr"), "dispIncr")

    self._info.log("Updating integrators states.")
    for integrator in self.integrators:
      integrator.updateState(t, disp)

    self._info.log("Outputting elastic solution.")
    for output in self.output.bin:
      output.writeField(t, self._istep, disp, self.solnField['label'])
    self._istep += 1      
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def pde_formulation():
  """
  Factory associated with Implicit.
  """
  return Implicit()


# End of file 
