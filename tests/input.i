[Components]
[]

[AuxVariables]
  [T_cell]
    order = FIRST
    family = LAGRANGE
    initial_condition = 10.0
  []
[]

[Executioner]
  type = Transient
  nl_rel_tol = 1
  l_tol = 1
  start_time = 0.0
  end_time = 2
  dt = 1.0
[]

[Outputs]
  exodus = true
[]
