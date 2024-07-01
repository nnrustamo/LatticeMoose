[Mesh]
    [gmg]
      type = GeneratedMeshGenerator # Can generate simple lines, rectangles and rectangular prisms
      dim = 2 # Dimension of the mesh
      nx = 10 # Number of elements in the x direction
      ny = 10 # Number of elements in the y direction
      xmax = 9.9 # Length of test chamber
      ymax = 9.9 # Test chamber radius
    []
    allow_renumbering = false
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
    execute_on = 'timestep_begin'
    nx = 10
    ny = 10
    q = 9
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
  num_steps = 400
[]

[Outputs]
    exodus = true # Output Exodus format
[]
