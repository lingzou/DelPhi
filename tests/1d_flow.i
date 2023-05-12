[FluidProperties]
  [eos]
    type = HeliumFluidProperties
  []
[]

[Components]
  [inlet]
    type = TDJ
    input = 'pipe(in)'
    v_bc = 1.0
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
[]

[Outputs]
  perf_graph = true
  exodus = true
  file_base = 1d_flow_flibe
[]
