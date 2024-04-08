[Mesh]
  type = FileMesh
  file = untitled.msh
[]

[FluidProperties]
  [co2]
    type = CO2FluidProperties
  []
[]

[UserObjects]
  [./eos]
    type = MoskitoEOS1P_FPModule
    fp = co2
  [../]    
  #[./viscosity]
  #  type = MoskitoViscosityWaterVogel
  #[../]
  [./viscosity]
    type = MoskitoEOS1P_FPVIS
    fp = co2
  [../]
[]

[Materials]
  [./left_vertical]
    type = MoskitoFluidWell_1p1c
    pressure = p
    temperature = T
    flowrate = q
    well_direction = -z
    well_type = 1
    eos_uo = eos
    viscosity_uo = viscosity
    well_diameter = 0.2032
    roughness_type = rough
    roughness = 1e-04
    gravity = '0 0 -9.8'
    block = left_vertical
    outputs = exodus
  [../]
  [./lat-left_vertical]
    type = MoskitoLatHeat_Inc_Formation_1p
    temperature_inner = T
    outer_diameters = '0.2032 0.219075 0.2286'
    conductivities = '100 0.7'
    formation_density = 2400
    formation_heat_capacity = 1000
    formation_temperature_function = grad_func
    convective_thermal_resistance = true
    nondimensional_time_function = Hasan_Kabir_2012
    outputs = exodus
    block = left_vertical
    formation_thermal_conductivity = 2.5
  [../]
  [./right_vertical_up]
    type = MoskitoFluidWell_1p1c
    pressure = p
    temperature = T
    flowrate = q
    well_direction = -z
    well_type = -1
    eos_uo = eos
    viscosity_uo = viscosity
    well_diameter = 0.2032
    roughness_type = rough
    roughness = 1e-04
    gravity = '0 0 -9.8'
    block = right_vertical_up
    outputs = exodus
  [../]
  [./lat_right_vertical_up]
    type = MoskitoLatHeat_Inc_Formation_1p
    temperature_inner = T
    outer_diameters = '0.2032 0.219075 0.2286'
    conductivities = '100 0.02'
    formation_density = 2400
    formation_heat_capacity = 1000
    formation_temperature_function = grad_func
    convective_thermal_resistance = true
    nondimensional_time_function = Hasan_Kabir_2012
    outputs = exodus
    block = right_vertical_up
    formation_thermal_conductivity = 2.5
  [../]
  [./right_vertical_bet]
    type = MoskitoFluidWell_1p1c
    pressure = p
    temperature = T
    flowrate = q
    well_direction = -z
    well_type = -1
    eos_uo = eos
    viscosity_uo = viscosity
    well_diameter = 0.2032
    roughness_type = rough
    roughness = 1e-04
    gravity = '0 0 -9.8'
    block = right_vertical_bet
    outputs = exodus
  [../]
  [./lat_right_vertical_bet]
    type = MoskitoLatHeat_Inc_Formation_1p
    temperature_inner = T
    outer_diameters = '0.2032 0.219075 0.2286'
    conductivities = '100 0.02'
    formation_density = 2400
    formation_heat_capacity = 1000
    formation_temperature_function = grad_func
    convective_thermal_resistance = true
    nondimensional_time_function = Hasan_Kabir_2012
    outputs = exodus
    block = right_vertical_bet
    formation_thermal_conductivity = 2.5
  [../]
  [./right_vertical_down]
    type = MoskitoFluidWell_1p1c
    pressure = p
    temperature = T
    flowrate = q
    well_direction = -z
    well_type = -1
    eos_uo = eos
    viscosity_uo = viscosity
    well_diameter = 0.2032
    roughness_type = rough
    roughness = 1e-04
    gravity = '0 0 -9.8'
    block = right_vertical_down
    outputs = exodus
  [../]
  [./lat_right_vertical_down]
    type = MoskitoLatHeat_Inc_Formation_1p
    temperature_inner = T
    outer_diameters = '0.2032 0.219075 0.2286'
    conductivities = '100 0.7'
    formation_density = 2400
    formation_heat_capacity = 1000
    formation_temperature_function = grad_func
    convective_thermal_resistance = true
    nondimensional_time_function = Hasan_Kabir_2012
    outputs = exodus
    block = right_vertical_down
    formation_thermal_conductivity = 2.5
  [../]
  [./horizontal]
     type = MoskitoFluidWell_1p1c
     pressure = p
     temperature = T
     flowrate = q
     well_direction = x
     well_type = 1
     eos_uo = eos
     viscosity_uo = viscosity
     well_diameter = 0.2032
     roughness_type = rough
     roughness = 2e-04
     gravity = '0 0 0'
     block = horizontal
     outputs = exodus
   [../]
   [./lat-horizontal]
     type = MoskitoLatHeat_Inc_Formation_1p
     temperature_inner = T
     outer_diameters = '0.2032 0.203201'
     conductivities = '3.000'
     formation_density = 2400
     formation_heat_capacity = 1000
     formation_temperature_function = grad_func
     convective_thermal_resistance = true
     nondimensional_time_function = Hasan_Kabir_2012
     outputs = exodus
     block = horizontal
     formation_thermal_conductivity = 2.5
   [../]
[]

[Functions]
  [./dts]
    type = PiecewiseLinear
    x = '0 7200  72000.0 943000  315360000 3153600000'
    y = '1000 1000.0 1000.0 360000 3600000 3600000'
  [../]
  [./conc]
    type = ParsedFunction
    value = 'if(t>16000 & t<73000, 1.72, 0.25)'
    # value = 0.4
  [../]
  [./grad_func]
    type = ParsedFunction
    value ='283.15 - z * 0.03'
  [../]
  [./conductivity_gradient]
    type = ParsedFunction
    value ='1.5 - z * 0.00075'
  [../]
[]

[BCs]
  [./pleft]
    type = DirichletBC
    variable = p
    boundary = inlet
    value = 10000000
  [../]
  [./qleft]
    type = DirichletBC
    variable = q
    boundary = inlet
    value = 0.02
  [../]
  [./tleft]
    type = DirichletBC
    variable = T
    boundary = inlet
    value = 303.15
  [../]
[]

[Variables]
  [./T]
    [./InitialCondition]
      type = FunctionIC
      function = '303.15 - z * 0.03'
      variable = T
    [../]
  [../]
  [./p]
    [./InitialCondition]
      type = FunctionIC
      function = '-10000 * z + 10000000'
      variable = p
    [../]
  [../]
  [./q]
    initial_condition = 0.02
  [../]
[]

[Kernels]
  [./Tkernel]
    type = MoskitoEnergy_1p1c
    variable = T
    flowrate = q
    pressure = p
    block = 'right_vertical_up right_vertical_bet right_vertical_down horizontal left_vertical'
    gravity_energy = false
  [../]
  [./Ttkernel]
    type = MoskitoTimeEnergy_1p1c
    variable = T
    flowrate = q
    pressure = p
    block = 'right_vertical_up right_vertical_bet right_vertical_down horizontal left_vertical'
  [../]
  [./pkernel]
    type = MoskitoMass_1p1c
    variable = p
    flowrate = q
    temperature = T
    block = 'right_vertical_up right_vertical_bet right_vertical_down horizontal left_vertical'
  [../]
  [./ptkernel]
    type = MoskitoTimeMass_1p1c
    variable = p
    temperature = T
    block = 'right_vertical_up right_vertical_bet right_vertical_down horizontal left_vertical'
  [../]
  [./qkernel]
    type = MoskitoMomentum_1p1c
    variable = q
    pressure = p
    temperature = T
    block = 'right_vertical_up right_vertical_bet right_vertical_down horizontal left_vertical'
  [../]
  [./qtkernel]
    type = MoskitoTimeMomentum_1p1c
    variable = q
    pressure = p
    temperature = T
    block = 'right_vertical_up right_vertical_bet right_vertical_down horizontal left_vertical'
  [../]
  [./heat]
    type = MoskitoLatHeatIncFormation_1p
    variable = T
    block = 'right_vertical_up right_vertical_bet right_vertical_down horizontal left_vertical'
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
  end_time = 315360
  l_max_its = 100
  solve_type = NEWTON
  steady_state_detection = true
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
  file_base = out_55
  exodus = true
[]
