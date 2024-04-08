[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 150
    ny = 100
    xmin = -2500
    xmax = 3500
    ymin = -1300
    ymax = 1300
  []
  [injection_node]
    input = gen
    type = ExtraNodesetGenerator
    new_boundary = injection_node
    coord = '-500 0 0'
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]

[Variables]
  [p]
    initial_condition = 50E6
  []
  [T]
    initial_condition = 433.15
    scaling = 1E-6 # fluid enthalpy is roughly 1E6
  []
[]

 [BCs]
  #  [temperature]
  #    type = DirichletBC
  #    variable = temperature
  #    value = 423.15
  #    boundary = 'left right top bottom'
  #  []
  #  [pressure]
  #   type = FunctionDirichletBC
  #   variable = porepressure
  #   function = csv_data
  #   boundary = injection_node
  # []
  # [ttt]
  #   type = FunctionDirichletBC
  #   variable = temperature
  #   function = csv_datat
  #   boundary = injection_node
  # []
  [./ttt]
    type = PostprocessorDirichletBC
    variable = T
    boundary = injection_node
    postprocessor = T_receiver2
  [../]
    [./ppp]
      type = PostprocessorDirichletBC
      variable = p
      boundary = injection_node
      postprocessor = P_receiver2
    [../]
 []

[Functions]
  [mass_flux_in_func]
    type = ParsedFunction
    expression = '10'
  []
  [mass_flux_out_func]
    type = ParsedFunction
    expression = '-10'
  []
  [T_in_func]
    type = ParsedFunction
    expression = '313.15'
  []
  [./csv_data]
    type = Moskito2Constant
    data_file = pcres.csv
    direction = RIGHT
    format = columns
    execute_on = TIMESTEP_BEGIN
  [../]
  [./csv_datat]
    type = Moskito2Constant
    data_file = tcres.csv
    direction = RIGHT_INCLUSIVE
    format = columns
    execute_on = TIMESTEP_BEGIN
  [../]
  [./csv_dataq]
    type = Moskito2Constant
    data_file = qcres.csv
    direction = RIGHT_INCLUSIVE
    format = columns
    execute_on = TIMESTEP_BEGIN
  [../]
  [./dts]
    type = PiecewiseLinear
    x = '0 10000  315360  31536000 630720000'
    y = '1000.0 1000.0 86400 3153600 1000000'
    #x = '0 7200  72000.0 943000  315360000 3153600000'
    #y = '1000 1000.0 1000.0 360000 3600000 3600000'
  [../]
[]

[Kernels]
  [mass_dot]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = p
  []
  [mass_flux]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = p
  []
  [energy_dot]
    type = PorousFlowEnergyTimeDerivative
    variable = T
  []
  [heat_advection]
    type = PorousFlowHeatAdvection
    variable = T
  []
  [heat_conduction]
    type = PorousFlowHeatConduction
    variable = T
  []
[]

[AuxVariables]
  [./qq_aux]
  [../]
    # [./tt_aux]
    #   [../]
[]

[AuxKernels]
  [./qq_aux2]
    type = ConstantAux
    variable = qq_aux
    value = 0.01
  [../]
    # [./pp_aux2]
    #   type = ConstantAux
    #   variable = pp_aux
    #   value = 30000000
    # [../]
[]

[Postprocessors]
  [mass_flux_in]
    type = FunctionValuePostprocessor
    function = mass_flux_in_func
    execute_on = 'initial timestep_begin'
  []
  [mass_flux_out]
    type = FunctionValuePostprocessor
    function = mass_flux_out_func
    execute_on = 'initial timestep_begin'
  []
  [T_in]
    type = FunctionValuePostprocessor
    function = T_in_func
    execute_on = 'initial timestep_begin'
  []
  [./temp_pro]
    type = PointValue
    point = '1500 0 0'
    variable = 'T'
  [../]
  [belphip]
      type = PointValue
      variable = p
      point = '1500 0 0'
  []
  [belphit]
      type = PointValue
      variable = T
      point = '1500 0 0'
  []

  [./P_value2]
    type = PointValue
    variable = p
    point = '1500 0 0'
  [../]
  [./T_value2]
    type = PointValue
    variable = T
    point = '1500 0 0'
  [../]
  [./Q_value2]
    type = PointValue
    variable = qq_aux
    point = '1500 0 0'
  [../]
  [./P_receiver2]
    type = Receiver
    # execute_on = NONLINEAR
  [../]
  [./T_receiver2]
    type = Receiver
    # execute_on = NONLINEAR
  [../]
  [./Q_receiver2]
    type = Receiver
    # execute_on = NONLINEAR
  [../]  
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'T p'  
    number_fluid_phases = 1
    number_fluid_components = 1
  []
  [pc]
    type = PorousFlowCapillaryPressureConst
  []
	[injected_mass]
    type = PorousFlowSumQuantity
  []
  [produced_mass]
    type = PorousFlowSumQuantity
  []
  [produced_heat]
    type = PorousFlowSumQuantity
  []
    [produced_mass_H2O1]
    type = PorousFlowSumQuantity
  []

  [produced_heat1]
    type = PorousFlowSumQuantity
  []
[]

[FluidProperties]
  [the_simple_fluid]
    type = SimpleFluidProperties
    thermal_expansion = 7E-4
    bulk_modulus = 2.5E9
    viscosity = 5E-4
    density0 = 1000
    cv = 4000.0
    cp = 4000.0
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = T
  []
  [phase]
    type = PorousFlow1PhaseP 
    porepressure = p
    capillary_pressure = pc
  []
  [massfrac]
    type = PorousFlowMassFraction
  []
  [water]
    type = PorousFlowSingleComponentFluid
    fp = the_simple_fluid
    phase = 0
  []
  [porosity]
    type=PorousFlowPorosityConst
    porosity = 0.1  
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1e-12 0 0  0 1e-12 0  0 0 1e-12' 
  []
  [relperm]
    type = PorousFlowRelativePermeabilityCorey
    n = 0
    phase = 0
  []
  [thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '2.8 0 0  0 2.8 0  0 0 2.8' 
  []
  [internal_energy]
    type = PorousFlowMatrixInternalEnergy
    specific_heat_capacity = 1000 
    density = 2500 
  []
  [undrained_density]
    type = GenericConstantMaterial
    prop_names = density
    prop_values = 2500 
  [] 
[]

[Preconditioning]
  active = basic
  [basic]
    type = SMP
    full = true
    petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix'
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = ' asm      lu           NONZERO                   2'
  []
  [preferred_but_might_not_be_installed]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       mumps'
  []
[]

[Executioner]
  type = Transient
  start_time = 0
  end_time = 1576800
  l_max_its = 100
  nl_abs_tol = 1e-6
  solve_type = NEWTON
  steady_state_detection = false
  [./TimeStepper]
    type = FunctionDT
    function = dts
  [../]
[]

[Outputs]
  file_base = pf_out
  exodus = true
  csv = true
[]

[DiracKernels]
  [inject_mass]
    type = PorousFlowPointSourceFromPostprocessor
    variable = p
    mass_flux = mass_flux_in
    point = '-500 0 0'
  []
  # [inject_heat]
  #   type = PorousFlowPointEnthalpySourceFromPostprocessor
  #   variable = T
  #   mass_flux = mass_flux_in
  #   point = '-500 0 0'
  #   T_in = T_in
  #   pressure = porepressure
  #   fp = the_simple_fluid
  # []

  [produce_H2O_1]
    type = PorousFlowPolyLineSink
    SumQuantityUO = produced_mass_H2O1
    fluxes = 10
    p_or_t_vals = 0.0
    line_length = 1.0
    point_file = /Users/morteza/Documents/projects/moose/modules/porous_flow/1.bh
    variable = p
  []

  [remove_heat_at_production_well_1]
    type = PorousFlowPolyLineSink
    SumQuantityUO = produced_heat1
    fluxes = 10
    p_or_t_vals = 0.0
    line_length = 1.0
    use_enthalpy = true
    point_file = /Users/morteza/Documents/projects/moose/modules/porous_flow/1.bh
    variable = T
  []
[]


# [Executioner]
#   type = Transient
#   solve_type = Newton
#   end_time = 900000000
#   # dt = 2E5
#   dt = 2E5
#   nl_abs_tol = 1e-9
# []


# [Postprocessors]
#   [belphip]
#     type = PointValue
#     variable = porepressure
#     point = '1500 0 0'
#     outputs = 'csv'
#   []
#   [belphit]
#     type = PointValue
#     variable = temperature
#     point = '1500 0 0'
#     outputs = 'csv'
#   []
# []

# [Postprocessors]
#   [belphip]
#     type = PointValue
#     variable = porepressure
#     point = '1500 0 0'
#     outputs = 'csv'
#   []
#   [belphit]
#     type = PointValue
#     variable = temperature
#     point = '1500 0 0'
#     outputs = 'csv'
#   []
# []

[MultiApps]
  [./coupled2]
    positions = '1500 0 0'
    type = TransientMultiApp
    input_files = /Users/morteza/Documents/projects/moskito2/production.i
    app_type = MoskitoApp
    execute_on = TIMESTEP_END
    sub_cycling = true
  [../]
[]

[Transfers]
  [P_transfer2]
    type = MultiAppPostprocessorTransfer
    to_multi_app = coupled2
    from_postprocessor = P_value2
    to_postprocessor = P_receiver3
  []
  [T_transfer2]
    type = MultiAppPostprocessorTransfer
    to_multi_app = coupled2
    from_postprocessor = T_value2
    to_postprocessor = T_receiver3
  []
  [Q_transfer2]
    type = MultiAppPostprocessorTransfer
    to_multi_app = coupled2
    from_postprocessor = Q_value2
    to_postprocessor = Q_receiver3
  []
[]