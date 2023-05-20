## Put two very different fluids in the same system, testing the numerical robustness

[FluidProperties]
  [helium]
    type = HeliumFluidProperties
  []
  [flibe]
    type = FlibeFluidProperties
  []
[]

[Components]
  [inlet-1]
    type = TDJ
    input = 'pipe-1(in)'
    v_bc = 1.0
    T_bc = 300.0
    eos = helium
  []

  [pipe-1]
    type = TestOneDFlow
    position = '0 0 0'
    orientation = '0 0 1'
    n_elems = 10
    length = 0.5
    eos = helium
    initial_P = 1e5
    initial_V = 0.5
    initial_T = 300
    Dh = 0.01
    f = 0.1
    A = 1
    heat_source = 1e5
  []

  [outlet-1]
    type = TDV
    input = 'pipe-1(out)'
    p_bc = 1.0e5
    T_bc = 300.0
    eos = helium
  []

  [inlet-2]
    type = TDJ
    input = 'pipe-2(in)'
    v_bc = 1.0
    T_bc = 600.0
    eos = flibe
  []

  [pipe-2]
    type = TestOneDFlow
    position = '0.5 0 0'
    orientation = '0 0 1'
    n_elems = 10
    length = 0.6
    eos = flibe
    initial_P = 1e5
    initial_V = 0.5
    initial_T = 300
    Dh = 0.01
    f = 0.1
    A = 1
    heat_source = 1e5
  []

  [outlet-2]
    type = TDV
    input = 'pipe-2(out)'
    p_bc = 1.0e5
    T_bc = 600.0
    eos = flibe
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
[]

[Outputs]
  perf_graph = true
  exodus = true
[]
