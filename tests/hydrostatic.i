[FluidProperties]
  [eos]
    type = LeadFluidProperties
  []
[]

[Components]
  [inlet]
    type = TDJ
    input = 'pipe(in)'
    v_bc = 0.0
    T_bc = 300.0
    eos = eos
  []

  [pipe]
    type = TestOneDFlow
    position = '0 0 0'
    orientation = '0 0 1'
    n_elems = 20
    length = 0.5
    eos = eos
    initial_P = 1e5
    initial_V = 0.5
    initial_T = 300
    Dh = 0.01
    f = 0.1
    A = 1
    heat_source = 0
  []

  [outlet]
    type = TDV
    input = 'pipe(out)'
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
  end_time = 5
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
