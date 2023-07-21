[FluidProperties]
  [eos]
    type = LinearizedFluidProperties
    rho_0 = 2883.649
    p_0 = 1e5
    T_0 = 873.0
    drho_dp = 0
    drho_dT = -0.887
    cp = 1051
  []
[]

[Components]
  [pipe-bottom]
    type = OneDFlowChannel
    position = '0 0 0'
    orientation = '1 0 0'
    order = 2
    n_elems = 10
    length = 0.1
    eos = eos
    initial_P = 1e5
    initial_V = 1e-3
    initial_T = 873
    Dh = 1
    A = 1
    WF_user_option = Welander
    Welander_constant = 288.36
    Hw = 3e5
    HT_surface_area_density = 1.0
    Tw = 883
  []

  [snj-bottom-to-right]
    type = SingleJunction
    inputs = 'pipe-bottom(out)'
    outputs = 'pipe-right(in)'
    eos = eos
  []

  [pipe-right]
    type = OneDFlowChannel
    position = '0.1 0 0'
    orientation = '0 0 1'
    order = 2
    n_elems = 10
    length = 0.9
    eos = eos
    initial_P = 1e5
    initial_V = 1e-3
    initial_T = 873
    Dh = 1
    A = 1
    WF_user_option = Welander
    Welander_constant = 288.36
  []

  [brv]
    type = VolumeBranch
    # the connection inputs can be very flexible
    inputs = 'pipe-right(out)'
    outputs = 'pipe-top(in) pipe-tdv(in)'
    eos = eos
    initial_P = 1e5
    initial_T = 873
    volume = 1.e-3 #0.0
  []

  [pipe-top]
    type = OneDFlowChannel
    position = '0.1 0 0.9'
    orientation = '-1 0 0'
    order = 2
    n_elems = 10
    length = 0.1
    eos = eos
    initial_P = 1e5
    initial_V = 1e-3
    initial_T = 873
    Dh = 1
    A = 1
    WF_user_option = Welander
    Welander_constant = 288.36
    Hw = 3e5
    HT_surface_area_density = 1.0
    Tw = 863
  []

  [snj-top-to-left]
    type = SingleJunction
    inputs = 'pipe-top(out)'
    outputs = 'pipe-left(in)'
    eos = eos
  []

  [pipe-left]
    type = OneDFlowChannel
    position = '0 0 0.9'
    orientation = '0 0 -1'
    order = 2
    n_elems = 10
    length = 0.9
    eos = eos
    initial_P = 1e5
    initial_V = 1e-3
    initial_T = 873
    Dh = 1
    A = 1
    WF_user_option = Welander
    Welander_constant = 288.36
  []

  [snj-left-to-bottom]
    type = SingleJunction
    inputs = 'pipe-left(out)'
    outputs = 'pipe-bottom(in)'
    eos = eos
    monitor = true
  []

  [pipe-tdv]
    type = OneDFlowChannel
    position = '0.1 0 0.9'
    orientation = '1 0 0'
    order = 2
    n_elems = 5
    length = 0.4
    eos = eos
    initial_P = 1e5
    initial_V = 0
    initial_T = 873
    Dh = 1
    A = 1
    WF_user_option = Welander
    Welander_constant = 10000
  []

  [tdv]
    type = TDV
    input = 'pipe-tdv(out)'
    p_bc = 1.0e5
    T_bc = 873
    eos = eos
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  nl_rel_tol = 1
  l_tol = 1
  start_time = 0.0
  end_time = 100
  dt = 0.2
  dtmin = 0.001

  # LeadFluidProperties is completely incompressible, i.e., drho_dp = 0,
  # and thus, we would like PTESc understand we have zeros on diagonal
  # and PETSc should use this option
  petsc_options_iname = '-pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'nonzero 1e-8'
[]

[Outputs]
  perf_graph = true
  exodus = true
  csv = true
[]
