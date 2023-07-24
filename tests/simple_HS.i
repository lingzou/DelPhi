[SolidProperties]
  [const_solid]
    type = ThermalFunctionSolidProperties
    rho = 1
    cp = 1
    k = 1
  []
[]

[Components]
  [hs]
    type = SimpleHeatStructure
    position = '0 0 0'
    orientation = '1 0 0'
    length = 0.2
    width = 1
    elem_number_axial = 4
    elem_number_width = 20
    solid = const_solid

    # boundary condition
    T_bc_left = 0
    T_bc_right = 0
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  nl_rel_tol = 1
  l_tol = 1
  start_time = 0.0
  end_time = 5
  dt = 1.0
  dtmin = 0.001
[]

[Outputs]
  exodus = true
[]
