[Mesh]
  type = LBMesh 
  dim = 2 # Dimension of the mesh
  nx = 40 # Number of nodes in the x direction
  ny = 40 # Number of nodes in the y direction
  nz = 1 # Number of nodes in the z direction
  xmin = 0 # Minimum x-coordinate
  xmax = 39.9 # Maximum x-coordinate
  ymin = 0 # Minimum y-coordinate
  ymax = 39.9 # Maximum y-coordinate
  zmin = 0 # Minimum z-coordinate
  zmax = 0 # Maximum z-coordinate

[]

[Variables]
  [u]
  []
[]

[AuxVariables]
  [v]
    # family = MONOMIAL
    # order = CONSTANT
  []
[]

# set up LBM simulation
[UserObjects]
  [LBM]
    type = LatticeBoltzmann
    n_subcycles = 5
    tolerance = 1.0e-4
    initial_density = 1.0
    taus = 0.6
    execute_on = 'timestep_begin'
  []
[]

[AuxKernels]
  [v]
    type = LBMVelocityAux
    variable = v
    lbm_uo = LBM
    execute_on = 'timestep_end'
  []
[]

[Problem]
  kernel_coverage_check = false
[]

[Executioner]
  type = Transient
  num_steps = 1000
[]

[Outputs]
  exodus = true # Output Exodus format
[]
