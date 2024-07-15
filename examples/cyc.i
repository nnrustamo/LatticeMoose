# Flow around cylinder

[Mesh]
  type = LBMesh 
  dim = 2 # Dimension of the mesh
  nx = 800 # Number of nodes in the x direction
  ny = 128 # Number of nodes in the y direction
  xmin = 0
  xmax = 0.8
  ymin = 0
  ymax = 0.5
[]

[AuxVariables]
  [velocity_x]
    family = LAGRANGE
    order = FIRST
  []
  [velocity_y]
    family = LAGRANGE
    order = FIRST
  []
  [velocity_z]
    family = LAGRANGE
    order = FIRST
  []
  [rho]
    family = LAGRANGE
    order = FIRST
  []
[]

[Variables]
  [dummy]
  []
[]


# set up LBM simulation
[UserObjects]
  [LBM]
    type = LatticeBoltzmann
    n_subcycles = 10
    tolerance = 1.0e-10
    initial_density = 1.5
    inlet_density = 2.0
    outlet_density = 1.5
    taus = 1.0
    execute_on = 'timestep_begin'
    fBody = 0
  []
[]

[AuxKernels]
  [vx]
    type = LBMDataAux
    variable = velocity_x
    lbm_uo = LBM
    var_type = 'vel_x'
    execute_on = 'timestep_end'
  []
  [vy]
    type = LBMDataAux
    variable = velocity_y
    lbm_uo = LBM
    var_type = 'vel_y'
    execute_on = 'timestep_end'
  []
  [vz]
    type = LBMDataAux
    variable = velocity_z
    lbm_uo = LBM
    var_type = 'vel_z'
    execute_on = 'timestep_end'
  []
  [rho]
    type = LBMDataAux
    variable = rho
    lbm_uo = LBM
    var_type = 'density'
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
