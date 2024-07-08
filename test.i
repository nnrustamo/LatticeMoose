[Mesh]
    [gmg]
      type = GeneratedMeshGenerator # Can generate simple lines, rectangles and rectangular prisms
      dim = 2 # Dimension of the mesh
      nx = 100 # Number of elements in the x direction
      ny = 100 # Number of elements in the y direction
      xmax = 99.9 # Length of test chamber
      ymax = 99.9 # Test chamber radius
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


[UserObjects]
  [test]
      type = testuser
      nx = 100
      ny = 100
      execute_on = 'timestep_begin'
  []
[]



[Problem]
  kernel_coverage_check = false
[]

[Executioner]
  type = Transient
  num_steps = 1
[]

[Outputs]
    exodus = true # Output Exodus format
[]
