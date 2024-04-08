[Mesh]
  type = FileMesh
  file = injection.msh
[]

[FluidProperties]
  [co2]
    type = CO2FluidProperties
  []
[]

[UserObjects]
  [./eos]
    #type = MoskitoEOS1P_Brine
    type = MoskitoEOS1P_FPModule
    fp = co2
  [../]
  [./viscosity]
    #type = MoskitoViscosityWaterVogel
    type = MoskitoEOS1P_FPVIS
    fp = co2
  [../]
[]

[Materials]
  [./l0_430]
    type = MoskitoFluidWell_1p1c
    pressure = p
    temperature = T
    flowrate = q
    well_direction = z
    well_type = 1
    eos_uo = eos
    viscosity_uo = viscosity
    well_diameter = 0.1429
    roughness_type = rough
    roughness = 1e-04
    gravity = '0 0 -9.8'
    # block = l0_430
    outputs = exodus
  [../]
  [./lat-l0_430]
    type = MoskitoLatHeat_Inc_Formation_1p
    temperature_inner = T
    outer_diameters = '0.1429 0.1778 0.2168'
    conductivities = '43.75 0.7'
    formation_density = 2200
    formation_heat_capacity = 850
    formation_temperature_function = grad_func
    convective_thermal_resistance = true
    nondimensional_time_function = Hasan_Kabir_2012
    outputs = exodus
    # block = l0_430
    formation_thermal_conductivity = 3.8
  [../]
[]

[Functions]
  [./dts]
    type = PiecewiseLinear
    x = '0 10000  315360  31536000 630720000'
    y = '1000.0 1000.0 86400 3153600 1000000'
    #x = '0 7200  72000.0 943000  315360000 3153600000'
    #y = '1000 1000.0 1000.0 360000 3600000 3600000'
  [../]
  [./conc]
    type = ParsedFunction
    value = 'if(t>16000 & t<73000, 1.72, 0.25)'
    # value = 0.4
  [../]
  [./grad_func]
    type = ParsedFunction
    value ='300.0 - z * 0.029'
  [../]
  [./conductivity_gradient]
    type = ParsedFunction
    value ='1.5 - z * 0.00075'
  [../]
[]

[BCs]
  [./pleft]
    type = PostprocessorDirichletBC
    variable = p
    boundary = outlet
    postprocessor = P_receiver3
  [../]
    # [./pleft]
    #   type = DirichletBC
    #   variable = p
    #   boundary = outlet
    #   value = 34000000
    # [../]
  # [./qleft]
  #   type = FunctionDirichletBC
  #   variable = q
  #   boundary = outlet
  #   function = '0.008'
  # [../]
  [./qleft]
    type = PostprocessorDirichletBC
    variable = q
    boundary = outlet
    postprocessor = Q_receiver3
  [../]
  # [./tleft]
  #   type = DirichletBC
  #   variable = T
  #   boundary = outlet
  #   value = 305.0
  # [../]
  [./tleft]
    type = PostprocessorDirichletBC
    variable = T
    boundary = outlet
    postprocessor = T_receiver3
  [../]
[]

[Variables]
  [./T]
    [./InitialCondition]
      type = FunctionIC
      function = '300.0 - z * 0.029'
      variable = T
    [../]
  [../]
  [./p]
    [./InitialCondition]
      type = FunctionIC
      function = '24000000 -10000 * z'
      variable = p
    [../]
  [../]
  [./q]
    initial_condition = 0.008
  [../]
[]

[Kernels]
  [./Tkernel]
    type = MoskitoEnergy_1p1c
    variable = T
    flowrate = q
    pressure = p
    #block = 'right_vertical_up right_vertical_bet right_vertical_down horizontal left_vertical'
    gravity_energy = false
  [../]
  [./Ttkernel]
    type = MoskitoTimeEnergy_1p1c
    variable = T
    flowrate = q
    pressure = p
    #block = 'right_vertical_up right_vertical_bet right_vertical_down horizontal left_vertical'
  [../]
  [./pkernel]
    type = MoskitoMass_1p1c
    variable = p
    flowrate = q
    temperature = T
    #block = 'right_vertical_up right_vertical_bet right_vertical_down horizontal left_vertical'
  [../]
  [./ptkernel]
    type = MoskitoTimeMass_1p1c
    variable = p
    temperature = T
    #block = 'right_vertical_up right_vertical_bet right_vertical_down horizontal left_vertical'
  [../]
  [./qkernel]
    type = MoskitoMomentum_1p1c
    variable = q
    pressure = p
    temperature = T
    #block = 'right_vertical_up right_vertical_bet right_vertical_down horizontal left_vertical'
  [../]
  [./qtkernel]
    type = MoskitoTimeMomentum_1p1c
    variable = q
    pressure = p
    temperature = T
    #block = 'right_vertical_up right_vertical_bet right_vertical_down horizontal left_vertical'
  [../]
  [./heat]
    type = MoskitoLatHeatIncFormation_1p
    variable = T
    #block = 'right_vertical_up right_vertical_bet right_vertical_down horizontal left_vertical'
  [../]
[]

[AuxVariables]
  [./HeatFlux]
      family = MONOMIAL
  [../]
[]


[Preconditioning]
  active = pn1
  [./p1]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
    petsc_options_value = ' bjacobi  ilu          NONZERO                 '
  [../]
  [./pn1]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -snes_type -snes_linesearch_type'
    petsc_options_value = ' bjacobi  ilu          NONZERO                   newtonls   basic               '
  [../]
[]

[Executioner]
  type = Transient
  start_time = 0
  end_time = 630720000
  l_max_its = 100
  solve_type = NEWTON
  steady_state_detection = false
  nl_abs_tol = 1e-6
  [./TimeStepper]
    type = FunctionDT
    function = dts
  [../]
[]

# [Outputs]
#   exodus = true
#   print_linear_residuals = true
# []

[Outputs]
  # execute_on = 'timestep_end'
  file_base = production_out
  exodus = true
[]

[Postprocessors]
  [./P_receiver3]
    type = Receiver
    # execute_on = NONLINEAR
  [../]
  [./T_receiver3]
    type = Receiver
    # execute_on = NONLINEAR
  [../]
  [./Q_receiver3]
    type = Receiver
    # execute_on = NONLINEAR
  [../]   
[]
