[Mesh]
  type = LBMesh 
  dim = 3 # Dimension of the mesh
  nx = 40 # Number of nodes in the x direction
  ny = 40 # Number of nodes in the y direction
  nz = 40 # Number of nodes in the z direction
  xmin = 0 # Minimum x-coordinate
  xmax = 39.9 # Maximum x-coordinate
  ymin = 0 # Minimum y-coordinate
  ymax = 39.9 # Maximum y-coordinate
  zmin = 0 # Minimum z-coordinate
  zmax = 39.9# Maximum z-coordinate
[]

[Variables]
  [u]
  []
[]

[AuxVariables]
  [v]
    family = MONOMIAL
    order = CONSTANT
  []
  [rho]
    family = MONOMIAL
    order = CONSTANT
  []
[]

# set up LBM simulation
[UserObjects]
  [LBM]
    type = LatticeBoltzmann
    n_subcycles = 50
    tolerance = 1.0e-4
    initial_density = 1.0
    inlet_density = 1.0
    outlet_density = 0.9
    taus = 0.6
    execute_on = 'timestep_begin'
    fBody = 0.0
  []
[]

[AuxKernels]
  [v]
    type = LBMVelocityAux
    variable = v
    lbm_uo = LBM
    execute_on = 'timestep_end'
  []
  [rho]
    type = LBMPressureAux
    variable = rho
    lbm_uo = LBM
    execute_on = 'timestep_end'
  []
[]


[Problem]
  kernel_coverage_check = false
[]

[Executioner]
  type = Transient
  num_steps = 50
[]

[Outputs]
  exodus = true # Output Exodus format
[]
