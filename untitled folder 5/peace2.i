[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 50
    ny = 50
    nz = 1
    xmin = -50
    xmax = 50
    ymin = -50
    ymax = 50
    zmin = 0
    zmax = 10
  []
  [central_nodes]
    input = gen
    type = ExtraNodesetGenerator
    new_boundary = central_nodes
    coord = '0 0 0; 0 0 10'
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Variables]
  [porepressure]
    initial_condition = 20E6
  []
  [temperature]
    initial_condition = 323
    scaling = 1E-6 # fluid enthalpy is roughly 1E6
  []
[]

[BCs]
  # [injection_temperature]
  #   type = DirichletBC
  #   variable = temperature
  #   value = 305
  #   boundary = central_nodes
  # []
  # [./all]
  #   type = FunctionDirichletBC
  #   variable = temperature
  #   boundary = central_nodes
  #   function = csv_data
  # [../]
    [./all]
      type = PostprocessorDirichletBC
      variable = temperature
      boundary = central_nodes
      postprocessor = T_receiver2
    [../]

[]

[DiracKernels]
  [fluid_injection]
    type = PorousFlowPeacemanBorehole
    variable = porepressure
    SumQuantityUO = injected_mass
    point_file = /Users/morteza/Documents/projects/moskito2/injection.bh
    function_of = pressure
    fluid_phase = 0
    bottom_p_or_t = 21E6
    unit_weight = '0 0 0'
    use_mobility = true
    character = -1
  []
  [fluid_production]
    type = PorousFlowPeacemanBorehole
    variable = porepressure
    SumQuantityUO = produced_mass
    point_file = /Users/morteza/Documents/projects/moskito2/production.bh
    function_of = pressure
    fluid_phase = 0
    bottom_p_or_t = 20E6
    unit_weight = '0 0 0'
    use_mobility = true
    character = 1
  []
  [remove_heat_at_production_well]
    type = PorousFlowPeacemanBorehole
    variable = temperature
    SumQuantityUO = produced_heat
    point_file = /Users/morteza/Documents/projects/moskito2/production.bh
    function_of = pressure
    fluid_phase = 0
    bottom_p_or_t = 20E6
    unit_weight = '0 0 0'
    use_mobility = true
    use_enthalpy = true
    character = 1
  []
[]

[UserObjects]
  [injected_mass]
    type = PorousFlowSumQuantity
  []
  [produced_mass]
    type = PorousFlowSumQuantity
  []
  [produced_heat]
    type = PorousFlowSumQuantity
  []
[]

[Postprocessors]
  [heat_joules_extracted_this_timestep]
    type = PorousFlowPlotQuantity
    uo = produced_heat
  []
[]

[FluidProperties]
  [the_simple_fluid]
    type = SimpleFluidProperties
    thermal_expansion = 2E-4
    bulk_modulus = 2E9
    viscosity = 1E-3
    density0 = 1000
    cv = 4000.0
    cp = 4000.0
  []
[]

[PorousFlowUnsaturated]
  porepressure = porepressure
  temperature = temperature
  coupling_type = ThermoHydro
  gravity = '0 0 0'
  fp = the_simple_fluid
[]

[Materials]
  [porosity]
    type = PorousFlowPorosityConst # only the initial value of this is ever used
    porosity = 0.1
  []
  [biot_modulus]
    type = PorousFlowConstantBiotModulus
    solid_bulk_compliance = 1E-10
    fluid_bulk_modulus = 2E9
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1E-12 0 0   0 1E-12 0   0 0 1E-12'
  []
  [thermal_expansion]
    type = PorousFlowConstantThermalExpansionCoefficient
    fluid_coefficient = 5E-6
    drained_coefficient = 2E-4
  []
  [thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '1 0 0  0 1 0  0 0 1'
  []
  [rock_heat]
    type = PorousFlowMatrixInternalEnergy
    density = 2500.0
    specific_heat_capacity = 1200.0
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

[Functions]
  [./dts]
    type = PiecewiseLinear
    x = '0 10000  315360  31536000 630720000'
    y = '1000.0 1000.0 86400 3153600 1000000'
    #x = '0 7200  72000.0 943000  315360000 3153600000'
    #y = '1000 1000.0 1000.0 360000 3600000 3600000'
  [../]
  # [./csv_data]
  #   type = MoskitoConstant
  #   data_file = pc.csv
  #   direction = RIGHT
  #   format = columns
  #   execute_on = TIMESTEP_BEGIN
  # [../]   
[]

[Executioner]
  type = Transient
  start_time = 0
  end_time = 630720000
  l_max_its = 100
  nl_abs_tol = 1e-11
  solve_type = NEWTON
  steady_state_detection = false
  [./TimeStepper]
    type = FunctionDT
    function = dts
  [../]
[]

# [Executioner]
#   type = Transient
#   solve_type = Newton
#   end_time = 2E6
#   dt = 2E5
# []

[Outputs]
  exodus = true
  file_base = peace_output
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
    value = 0.008
  [../]
    # [./pp_aux2]
    #   type = ConstantAux
    #   variable = pp_aux
    #   value = 30000000
    # [../]
[]

[Postprocessors]
  [./P_value2]
    type = PointValue
    variable = porepressure
    point = '0 0 10'
  [../]
  [./T_value2]
    type = PointValue
    variable = temperature
    point = '0 0 10'
  [../]
  [./Q_value2]
    type = PointValue
    variable = qq_aux
    point = '0 0 10'
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

[MultiApps]
  [./coupled2]
    positions = '0 0 10'
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