[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 30
    ny = 20
    nz = 20
    xmin = -1500
    xmax = 1500
    ymin = -1000
    ymax = 1000
    zmin = -100
    zmax = 0
  []
  [injection_node1]
    input = gen
    type = ExtraNodesetGenerator
    new_boundary = injection_node1
    coord = '-500 0 -30; -500 0 -35; -500 0 -40; -500 0 -45; -500 0 -50; -500 0 -55; -500 0 -60; -500 0 -65; -500 0 -70'
    # tolerance = 5
  []
  # [injection_node2]
  #   input = gen
  #   type = ExtraNodesetGenerator
  #   new_boundary = injection_node2
  #   coord = '-500 0 -35'
  # []
  # [injection_node3]
  #   input = gen
  #   type = ExtraNodesetGenerator
  #   new_boundary = injection_node3
  #   coord = '-500 0 -40'
  # []
  # [injection_node4]
  #   input = gen
  #   type = ExtraNodesetGenerator
  #   new_boundary = injection_node4
  #   coord = '-500 0 -45'
  # []
  # [injection_node5]
  #   input = gen
  #   type = ExtraNodesetGenerator
  #   new_boundary = injection_node5
  #   coord = '-500 0 -50'
  # []
  # [injection_node6]
  #   input = gen
  #   type = ExtraNodesetGenerator
  #   new_boundary = injection_node6
  #   coord = '-500 0 -55'
  # []
  # [injection_node7]
  #   input = gen
  #   type = ExtraNodesetGenerator
  #   new_boundary = injection_node7
  #   coord = '-500 0 -60'
  # []
  # [injection_node8]
  #   input = gen
  #   type = ExtraNodesetGenerator
  #   new_boundary = injection_node8
  #   coord = '-500 0 -65'
  # []
  # [injection_node9]
  #   input = gen
  #   type = ExtraNodesetGenerator
  #   new_boundary = injection_node9
  #   coord = '-500 0 -70'
  # []
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 -9.8'
[]

[Variables]
  [p]
    # initial_condition = 30E6
    [./InitialCondition]
      type = FunctionIC
      function = '50E6 - 1000*9.8*z'
      variable = T
    [../]
  []
  [T]
    # initial_condition = 383.15
    scaling = 1E-6 # fluid enthalpy is roughly 1E6
    [./InitialCondition]
      type = FunctionIC
      function = '293.15 + 5*30 - z*0.03'
      variable = T
    [../]
  []
[]

 [BCs]
   [T]
     type = DirichletBC
     variable = T
     value = 443.15
    #  boundary = 'left right front back top bottom'
     boundary = 'top'
   []
   [T2]
    type = DirichletBC
    variable = p
    value = 50E6
    boundary = 'top'
  []
  # [TTT]
  #   type = DirichletBC
  #   variable = T
  #   value = 333.15
  #   # boundary = 'injection_node1 injection_node2 injection_node3 injection_node4 injection_node5 injection_node6 injection_node7 injection_node8 injection_node9'
  #   boundary = 'injection_node1'
  # []
  # [ppp]
  #   type = DirichletBC
  #   variable = p
  #   value = 55000000
  #   # boundary = 'injection_node1 injection_node2 injection_node3 injection_node4 injection_node5 injection_node6 injection_node7 injection_node8 injection_node9'
  #   boundary = 'injection_node1'
  # []

  [./ttt]
    type = PostprocessorDirichletBC
    variable = T
    boundary = 'injection_node1'
    postprocessor = T_receiver2
  [../]
    [./ppp]
      type = PostprocessorDirichletBC
      variable = p
      boundary = 'injection_node1'
      postprocessor = P_receiver2
    [../]
  #  [pressure]
  #   type = DirichletBC
  #   variable = p
  #   value = 32000000
  #   boundary = injection_node
  # []
  # [pressure]
  #   type = DirichletBC
  #   variable = p
  #   value = 32000000
  #   boundary = 'left right front top bottom'
  # []
 []

[Functions]
  [mass_flux_in_func]
    type = ParsedFunction
    expression = '5.0'
  []
  [mass_flux_out_func]
    type = ParsedFunction
    expression = '-5.0'
  []
  [T_in_func]
    type = ParsedFunction
    expression = '333.15'
  []
  [./dts]
    type = PiecewiseLinear
    x = '0 10000  315360  31536000 630720000'
    y = '1000.0 1000.0 86400 3153600 1000000'
    #x = '0 7200  72000.0 943000  315360000 3153600000'
    #y = '1000 1000.0 1000.0 360000 3600000 3600000'
  [../]
[]


[DiracKernels]
  [inject_mass1]
    type = PorousFlowPointSourceFromPostprocessor
    variable = p
    mass_flux = mass_flux_in
    point = '-500 0 -30'
  []
  [inject_mass2]
    type = PorousFlowPointSourceFromPostprocessor
    variable = p
    mass_flux = mass_flux_in
    point = '-500 0 -35'
  []
  [inject_mass3]
    type = PorousFlowPointSourceFromPostprocessor
    variable = p
    mass_flux = mass_flux_in
    point = '-500 0 -40'
  []
  [inject_mass4]
    type = PorousFlowPointSourceFromPostprocessor
    variable = p
    mass_flux = mass_flux_in
    point = '-500 0 -45'
  []
  [inject_mass5]
    type = PorousFlowPointSourceFromPostprocessor
    variable = p
    mass_flux = mass_flux_in
    point = '-500 0 -50'
  []
  [inject_mass6]
    type = PorousFlowPointSourceFromPostprocessor
    variable = p
    mass_flux = mass_flux_in
    point = '-500 0 -55'
  []
  [inject_mass7]
    type = PorousFlowPointSourceFromPostprocessor
    variable = p
    mass_flux = mass_flux_in
    point = '-500 0 -60'
  []
  [inject_mass8]
    type = PorousFlowPointSourceFromPostprocessor
    variable = p
    mass_flux = mass_flux_in
    point = '-500 0 -65'
  []
  [inject_mass9]
    type = PorousFlowPointSourceFromPostprocessor
    variable = p
    mass_flux = mass_flux_in
    point = '-500 0 -70'
  []
  # [inject_heat1]
  #   type = PorousFlowPointEnthalpySourceFromPostprocessor
  #   variable = T
  #   mass_flux = mass_flux_in
  #   point = '-500 0 -30'
  #   T_in = T_in
  #   pressure = p
  #   fp = the_simple_fluid
  # []
  # [inject_heat2]
  #   type = PorousFlowPointEnthalpySourceFromPostprocessor
  #   variable = T
  #   mass_flux = mass_flux_in
  #   point = '-500 0 -35'
  #   T_in = T_in
  #   pressure = p
  #   fp = the_simple_fluid
  # []
  # [inject_heat3]
  #   type = PorousFlowPointEnthalpySourceFromPostprocessor
  #   variable = T
  #   mass_flux = mass_flux_in
  #   point = '-500 0 -40'
  #   T_in = T_in
  #   pressure = p
  #   fp = the_simple_fluid
  # []
  # [inject_heat4]
  #   type = PorousFlowPointEnthalpySourceFromPostprocessor
  #   variable = T
  #   mass_flux = mass_flux_in
  #   point = '-500 0 -45'
  #   T_in = T_in
  #   pressure = p
  #   fp = the_simple_fluid
  # []
  # [inject_heat5]
  #   type = PorousFlowPointEnthalpySourceFromPostprocessor
  #   variable = T
  #   mass_flux = mass_flux_in
  #   point = '-500 0 -50'
  #   T_in = T_in
  #   pressure = p
  #   fp = the_simple_fluid
  # []
  # [inject_heat6]
  #   type = PorousFlowPointEnthalpySourceFromPostprocessor
  #   variable = T
  #   mass_flux = mass_flux_in
  #   point = '-500 0 -55'
  #   T_in = T_in
  #   pressure = p
  #   fp = the_simple_fluid
  # []
  # [inject_heat7]
  #   type = PorousFlowPointEnthalpySourceFromPostprocessor
  #   variable = T
  #   mass_flux = mass_flux_in
  #   point = '-500 0 -60'
  #   T_in = T_in
  #   pressure = p
  #   fp = the_simple_fluid
  # []
  # [inject_heat8]
  #   type = PorousFlowPointEnthalpySourceFromPostprocessor
  #   variable = T
  #   mass_flux = mass_flux_in
  #   point = '-500 0 -65'
  #   T_in = T_in
  #   pressure = p
  #   fp = the_simple_fluid
  # []
  # [inject_heat9]
  #   type = PorousFlowPointEnthalpySourceFromPostprocessor
  #   variable = T
  #   mass_flux = mass_flux_in
  #   point = '-500 0 -70'
  #   T_in = T_in
  #   pressure = p
  #   fp = the_simple_fluid
  # []

  [produce_H2O_1]
    type = PorousFlowPolyLineSink
    SumQuantityUO = produced_mass_H2O1
    fluxes = 1
    p_or_t_vals = 0.0
    line_length = 5.0
    point_file = /Users/morteza/Documents/projects/moskito2/3.bh
    variable = p
  []

  [remove_heat_at_production_well_1]
    type = PorousFlowPolyLineSink
    SumQuantityUO = produced_heat1
    fluxes = 1
    p_or_t_vals = 0.0
    line_length = 5.0
    use_enthalpy = true
    point_file = /Users/morteza/Documents/projects/moskito2/3.bh
    variable = T
  []
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
    value = 0.045
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
    point = '500 0 0'
    variable = 'T'
  [../]
  [./P_value2]
    type = PointValue
    variable = p
    point = '500 0 -50'
  [../]
  [./T_value2]
    type = PointValue
    variable = T
    point = '500 0 -50'
  [../]
  [./Q_value2]
    type = PointValue
    variable = qq_aux
    point = '500 0 -50'
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
  [T]
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
    permeability = '1e-13 0 0  0 1e-13 0  0 0 1e-13' 
  []
  [relperm]
    type = PorousFlowRelativePermeabilityCorey
    n = 0
    phase = 0
  []
  [thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '2.5 0 0  0 2.5 0  0 0 2.5' 
  []
  [internal_energy]
    type = PorousFlowMatrixInternalEnergy
    specific_heat_capacity = 920 
    density = 2600 
  []
  [undrained_density]
    type = GenericConstantMaterial
    prop_names = density
    prop_values = 2400 
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

# [Executioner]
#   type = Transient
#   solve_type = Newton
#   end_time = 900000000
#   # dt = 2E5
#   dt = 2E5
#   nl_abs_tol = 1e-9
# []

# [Executioner]
#   type = Transient
#   [./TimeStepper]
#     type = IterationAdaptiveDT
#     dt = 3600
#     growth_factor = 2
#   [../]
#   dtmax = 864000
#   end_time = 935360000
#   solve_type = NEWTON
# []

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
  exodus = true
  file_base = reservoir
[]

[MultiApps]
  [./coupled2]
    positions = '500 0 -50'
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