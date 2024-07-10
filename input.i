[Mesh]
  type = LBMesh 
  dim = 2 # Dimension of the mesh
  nx = 500 # Number of nodes in the x direction
  ny = 500 # Number of nodes in the y direction
  load_mesh_from_file = true
  mesh_file = 'porous_media/porous_media.txt'
[]

[Variables]
  [u]
  []
[]

[AuxVariables]
  [v]
    family = LAGRANGE
    order = FIRST
  []
  [rho]
    family = LAGRANGE
    order = FIRST
  []
[]

# set up LBM simulation
[UserObjects]
  [LBM]
    type = LatticeBoltzmann
    n_subcycles = 10
    tolerance = 1.0e-10
    initial_density = 1.0
    inlet_density = 1.0
    outlet_density = 1.0
    taus = 0.6
    execute_on = 'timestep_begin'
    fBody = 1.0e-6
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
  num_steps = 10
[]

[Outputs]
  exodus = true # Output Exodus format
[]
