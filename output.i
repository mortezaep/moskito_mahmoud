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
  [./viscosity]
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
    well_direction = -z
    well_type = 1
    eos_uo = eos
    viscosity_uo = viscosity
    well_diameter = 0.1429
    roughness_type = rough
    roughness = 0.0001
    gravity = '0 0 -9.8'
    block = l0_430
    outputs = exodus
  [../]
  [./lat-l0_430]
    type = MoskitoLatHeat_Inc_Formation_1p
    temperature_inner = T
    outer_diameters = '0.1429 0.1778 0.2168 0.2450 0.3179 0.3397 0.4'
    conductivities = '43.75 0.7 43.75 0.7 43.75 0.7'
    formation_density = 2200.0
    formation_heat_capacity = 850.0
    formation_temperature_function = grad_func
    convective_thermal_resistance = true
    nondimensional_time_function = Hasan_Kabir_2012
    outputs = exodus
    block = l0_430
    formation_thermal_conductivity = 2.0
  [../]
  [./l430_2150]
    type = MoskitoFluidWell_1p1c
    pressure = p
    temperature = T
    flowrate = q
    well_direction = -z
    well_type = 1
    eos_uo = eos
    viscosity_uo = viscosity
    well_diameter = 0.1429
    roughness_type = rough
    roughness = 0.0001
    gravity = '0 0 -9.8'
    block = l430_2150
    outputs = exodus
  [../]
  [./lat_l430_2150]
    type = MoskitoLatHeat_Inc_Formation_1p
    temperature_inner = T
    outer_diameters = '0.1429 0.1778 0.2168 0.2450 0.3179'
    conductivities = '43.75 0.7 43.75 0.7'
    formation_density = 2200.0
    formation_heat_capacity = 850.0
    formation_temperature_function = grad_func
    convective_thermal_resistance = true
    nondimensional_time_function = Hasan_Kabir_2012
    outputs = exodus
    block = l430_2150
    formation_thermal_conductivity = 2.0
  [../]
  [./l2150_5000]
    type = MoskitoFluidWell_1p1c
    pressure = p
    temperature = T
    flowrate = q
    well_direction = -z
    well_type = 1
    eos_uo = eos
    viscosity_uo = viscosity
    well_diameter = 0.1429
    roughness_type = rough
    roughness = 0.0001
    gravity = '0 0 -9.8'
    block = l2150_5000
    outputs = exodus
  [../]
  [./lat_l2150_5000]
    type = MoskitoLatHeat_Inc_Formation_1p
    temperature_inner = T
    outer_diameters = '0.1429 0.1778 0.2168'
    conductivities = '43.75 0.7'
    formation_density = 2200.0
    formation_heat_capacity = 850.0
    formation_temperature_function = grad_func
    convective_thermal_resistance = true
    nondimensional_time_function = Hasan_Kabir_2012
    outputs = exodus
    block = l2150_5000
    formation_thermal_conductivity = 2.0
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
     well_diameter = 0.1429
     roughness_type = rough
     roughness = 2e-04
     gravity = '0 0 0'
     block = horizontal
     outputs = exodus
   [../]
   [./lat-horizontal]
     type = MoskitoLatHeat_Inc_Formation_1p
     temperature_inner = T
     outer_diameters = '0.1429 0.1529'
     conductivities = '2.0'
     formation_density = 2200.0
     formation_heat_capacity = 850.0
     formation_temperature_function = grad_func
     convective_thermal_resistance = true
     nondimensional_time_function = Hasan_Kabir_2012
     outputs = exodus
     block = horizontal
     formation_thermal_conductivity = 2.0
   [../]
   [./r2150_5000]
     type = MoskitoFluidWell_1p1c
     pressure = p
     temperature = T
     flowrate = q
     well_direction = -z
     well_type = -1
     eos_uo = eos
     viscosity_uo = viscosity
     well_diameter = 0.1429
     roughness_type = rough
     roughness = 0.0001
     gravity = '0 0 -9.8'
     block = r2150_5000
     outputs = exodus
    [../]
    [./lat_r2150_5000]
    type = MoskitoLatHeat_Inc_Formation_1p
     temperature_inner = T
     outer_diameters = '0.1429 0.1778 0.2168'
     conductivities = '43.75 0.7'
     formation_density = 2200.0
     formation_heat_capacity = 850.0
     formation_temperature_function = grad_func
     convective_thermal_resistance = true
     nondimensional_time_function = Hasan_Kabir_2012
     outputs = exodus
     block = r2150_5000
     formation_thermal_conductivity = 2.0
   [../]
   [./r430_2150]
     type = MoskitoFluidWell_1p1c
     pressure = p
     temperature = T
     flowrate = q
     well_direction = -z
     well_type = -1
     eos_uo = eos
     viscosity_uo = viscosity
     well_diameter = 0.1429
     roughness_type = rough
     roughness = 0.0001
     gravity = '0 0 -9.8'
     block = r430_2150
     outputs = exodus
   [../]
   [./lat_r430_2150]
     type = MoskitoLatHeat_Inc_Formation_1p
     temperature_inner = T
     outer_diameters = '0.1429 0.1778 0.2168 0.2450 0.3179'
     conductivities = '43.75 0.7 43.75 0.7'
     formation_density = 2200.0
     formation_heat_capacity = 850.0
     formation_temperature_function = grad_func
     convective_thermal_resistance = true
     nondimensional_time_function = Hasan_Kabir_2012
     outputs = exodus
     block = r430_2150
     formation_thermal_conductivity = 2.0
   [../]
    [./r0_430]
      type = MoskitoFluidWell_1p1c
      pressure = p
      temperature = T
      flowrate = q
      well_direction = -z
      well_type = -1
      eos_uo = eos
      viscosity_uo = viscosity
      well_diameter = 0.1429
      roughness_type = rough
      roughness = 0.0001
      gravity = '0 0 -9.8'
      block = r0_430
      outputs = exodus
    [../]
    [./lat_r0_430]
      type = MoskitoLatHeat_Inc_Formation_1p
      temperature_inner = T
      outer_diameters = '0.1429 0.1778 0.2168 0.2450 0.3179 0.3397 0.4'
      conductivities = '43.75 0.3 43.75 0.7 43.75 0.7'
      formation_density = 2200.0
      formation_heat_capacity = 850.0
      formation_temperature_function = grad_func
      convective_thermal_resistance = true
      nondimensional_time_function = Hasan_Kabir_2012
      outputs = exodus
      block = r0_430
      formation_thermal_conductivity = 2.0
    [../]
[]

[Functions]
  [./dts]
    type = PiecewiseLinear
    x = '0 10000  315360  31536000 630720000'
    y = '1000.0 1000.0 86400 3153600 1000000'
  [../]
  [./conc]
    type = ParsedFunction
    value = 'if(t>16000 & t<73000, 1.72, 0.25)'
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
    type = DirichletBC
    variable = p
    boundary = inlet
    value = 24000000
  [../]
  [./qleft]
    type = FunctionDirichletBC
    variable = q
    boundary = inlet
    function = '0.008'
  [../]
  [./tleft]
    type = DirichletBC
    variable = T
    boundary = inlet
    value = 280.0
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
    gravity_energy = false
  [../]
  [./Ttkernel]
    type = MoskitoTimeEnergy_1p1c
    variable = T
    flowrate = q
    pressure = p
  [../]
  [./pkernel]
    type = MoskitoMass_1p1c
    variable = p
    flowrate = q
    temperature = T
  [../]
  [./ptkernel]
    type = MoskitoTimeMass_1p1c
    variable = p
    temperature = T
  [../]
  [./qkernel]
    type = MoskitoMomentum_1p1c
    variable = q
    pressure = p
    temperature = T
  [../]
  [./qtkernel]
    type = MoskitoTimeMomentum_1p1c
    variable = q
    pressure = p
    temperature = T
  [../]
  [./heat]
    type = MoskitoLatHeatIncFormation_1p
    variable = T
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
  [./TimeStepper]
    type = FunctionDT
    function = dts
  [../]
[]

[Outputs]
  # execute_on = 'timestep_end'
  file_base = out_55
  exodus = true
[]