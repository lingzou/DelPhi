## Not the best physical setup, for testing numerical robustness
[GlobalSimParameters]
  gravity = '0 0 0'
[]

[FluidProperties]
  [eos]
    type = HeliumFluidProperties
  []
[]

[Components]
  [inlet]
    type = TDV
    input = 'pipe(in)'
    p_bc = 1.0e5
    T_bc = 300.0
    eos = eos
  []

  [pipe]
    type = OneDFlowChannel
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
    heat_source = 1e5
  []

  [outlet]
    type = TDJ
    input = 'pipe(out)'
    v_bc = 1.0
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
[]

[Outputs]
  perf_graph = true
  exodus = true
[]
