[FluidProperties]
  [eos]
    type = LeadFluidProperties
  []
[]

[Components]
  [inlet]
    type = TDJ
    input = 'pipe-1(in)'
    v_bc = 1.0
    T_bc = 650
    eos = eos
  []

  [pipe-1]
    type = TestOneDFlow
    position = '0 0 0'
    orientation = '0 0 1'
    n_elems = 10
    length = 0.5
    eos = eos
    initial_P = 1e5
    initial_V = 0.5
    initial_T = 650
    Dh = 0.01
    f = 0.1
    A = 1
    heat_source = 1e5
  []

  [brv]
    type = VolumeBranch
    # the connection inputs can be very flexible
    inputs = 'pipe-1(out) pipe-2(in) pipe-3(in)'
    outputs = ''
    eos = eos
    initial_P = 1e5
    initial_T = 650
    volume = 1.0
  []

  [pipe-2]
    type = TestOneDFlow
    position = '-0.1 0 0.5'
    orientation = '0 0 1'
    n_elems = 10
    length = 0.5
    eos = eos
    initial_P = 1e5
    initial_V = 0.5
    initial_T = 650
    Dh = 0.01
    f = 0.1
    A = 1
    heat_source = 1e5
  []

  [pipe-3]
    type = TestOneDFlow
    position = '0.1 0 0.5'
    orientation = '0 0 1'
    n_elems = 10
    length = 0.5
    eos = eos
    initial_P = 1e5
    initial_V = 0.5
    initial_T = 650
    Dh = 0.01
    f = 0.1
    A = 1
    heat_source = 1e5
  []

  [outlet-2]
    type = TDV
    input = 'pipe-2(out)'
    p_bc = 1.0e5
    T_bc = 300.0
    eos = eos
  []

  [outlet-3]
    type = TDV
    input = 'pipe-3(out)'
    p_bc = 1.0e5
    T_bc = 300.0
    eos = eos
  []
[]

[Executioner]
  type = Transient
  nl_rel_tol = 1
  l_tol = 1
  start_time = 0.0
  end_time = 10
  dt = 1
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
[]