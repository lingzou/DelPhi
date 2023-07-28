[GlobalSimParameters]
  gravity = '0 0 0'
[]

[FluidProperties]
  [eos]
    type = LeadFluidProperties
  []
[]

[SolidProperties]
  [const_solid]
    type = ThermalFunctionSolidProperties
    rho = 1000
    cp = 200
    k = 50
  []
[]

[Components]
  [hs]
    type = SimpleHeatStructure2
    position = '0 0 0.01'
    orientation = '1 0 0'
    length = 0.5
    width = 0.01
    elem_number_axial = 10
    elem_number_width = 5
    solid = const_solid
    Ts_init = 400

    # boundary condition
    # T_bc_left = 400
    name_comp_left = pipe-1
    HT_surface_area_density_left = 200
    Hw_left = 1000
    T_bc_right = 400
  []

  [inlet]
    type = TDJ
    input = 'pipe-1(in)'
    v_bc = 1.0
    T_bc = 650
    eos = eos
  []

  [pipe-1]
    type = OneDFlowChannel
    position = '0 0 0'
    orientation = '1 0 0'
    order = 2
    n_elems = 10
    length = 0.5
    eos = eos
    initial_P = 1e5
    initial_V = 0.0
    initial_T = 650
    Dh = 0.01
    f = 0.1
    A = 1e-4
    heat_source = 1e8
  []

  [outlet]
    type = TDV
    input = 'pipe-1(out)'
    p_bc = 1.0e5
    T_bc = 650.0
    eos = eos
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  nl_rel_tol = 1
  l_tol = 1
  start_time = 0.0
  end_time = 40
  dt = 2.0
  dtmin = 0.001

  petsc_options_iname = '-pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'nonzero 1e-8'
[]

[Outputs]
  perf_graph = true
  exodus = true
  file_base = HS_OneDFlow_CHT_option_1
[]
