[Components]
  [comp-1]
    type = TestComponent
    position = '0 0 0'
    orientation = '0 0 1'
    n_elems = 20
    length = 0.5
  []

  [comp-2]
    type = TestComponent
    position = '0 0 1'
    orientation = '1 0 0'
    n_elems = 10
    length = 1
  []
[]

[Executioner]
  type = Transient
  nl_rel_tol = 1
  l_tol = 1
  start_time = 0.0
  end_time = 10
  dt = 1.0
  dtmin = 0.001
[]

[Outputs]
  exodus = true
[]
