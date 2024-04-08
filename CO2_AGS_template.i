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
    well_diameter = {{ diameter }}
    roughness_type = rough
    roughness = {{ epsilon }}
    gravity = '0 0 -9.8'
    block = l0_430
    outputs = exodus
  [../]
  [./lat-l0_430]
    type = MoskitoLatHeat_Inc_Formation_1p
    temperature_inner = T
    outer_diameters = '{{ diameter }} 0.1778 0.2168 0.2450 0.3179 0.3397 0.4'
    conductivities = '43.75 0.7 43.75 0.7 43.75 0.7'
    formation_density = {{ rock_density }}
    formation_heat_capacity = {{ rock_cp }}
    formation_temperature_function = grad_func
    convective_thermal_resistance = true
    nondimensional_time_function = Hasan_Kabir_2012
    outputs = exodus
    block = l0_430
    formation_thermal_conductivity = {{ rock_k }}
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
    well_diameter = {{ diameter }}
    roughness_type = rough
    roughness = {{ epsilon }}
    gravity = '0 0 -9.8'
    block = l430_2150
    outputs = exodus
  [../]
  [./lat_l430_2150]
    type = MoskitoLatHeat_Inc_Formation_1p
    temperature_inner = T
    outer_diameters = '{{ diameter }} 0.1778 0.2168 0.2450 0.3179'
    conductivities = '43.75 0.7 43.75 0.7'
    formation_density = {{ rock_density }}
    formation_heat_capacity = {{ rock_cp }}
    formation_temperature_function = grad_func
    convective_thermal_resistance = true
    nondimensional_time_function = Hasan_Kabir_2012
    outputs = exodus
    block = l430_2150
    formation_thermal_conductivity = {{ rock_k }}
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
    well_diameter = {{ diameter }}
    roughness_type = rough
    roughness = {{ epsilon }}
    gravity = '0 0 -9.8'
    block = l2150_5000
    outputs = exodus
  [../]
  [./lat_l2150_5000]
    type = MoskitoLatHeat_Inc_Formation_1p
    temperature_inner = T
    outer_diameters = '{{ diameter }} 0.1778 0.2168'
    conductivities = '43.75 0.7'
    formation_density = {{ rock_density }}
    formation_heat_capacity = {{ rock_cp }}
    formation_temperature_function = grad_func
    convective_thermal_resistance = true
    nondimensional_time_function = Hasan_Kabir_2012
    outputs = exodus
    block = l2150_5000
    formation_thermal_conductivity = {{ rock_k }}
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
     well_diameter = {{ diameter }}
     roughness_type = rough
     roughness = 2e-04
     gravity = '0 0 0'
     block = horizontal
     outputs = exodus
   [../]
   [./lat-horizontal]
     type = MoskitoLatHeat_Inc_Formation_1p
     temperature_inner = T
     outer_diameters = '{{ diameter }} {{ diameter+0.01 }}'
     conductivities = '2.0'
     formation_density = {{ rock_density }}
     formation_heat_capacity = {{ rock_cp }}
     formation_temperature_function = grad_func
     convective_thermal_resistance = true
     nondimensional_time_function = Hasan_Kabir_2012
     outputs = exodus
     block = horizontal
     formation_thermal_conductivity = {{ rock_k }}
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
     well_diameter = {{ diameter }}
     roughness_type = rough
     roughness = {{ epsilon }}
     gravity = '0 0 -9.8'
     block = r2150_5000
     outputs = exodus
    [../]
    [./lat_r2150_5000]
    type = MoskitoLatHeat_Inc_Formation_1p
     temperature_inner = T
     outer_diameters = '{{ diameter }} 0.1778 0.2168'
     conductivities = '43.75 0.7'
     formation_density = {{ rock_density }}
     formation_heat_capacity = {{ rock_cp }}
     formation_temperature_function = grad_func
     convective_thermal_resistance = true
     nondimensional_time_function = Hasan_Kabir_2012
     outputs = exodus
     block = r2150_5000
     formation_thermal_conductivity = {{ rock_k }}
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
     well_diameter = {{ diameter }}
     roughness_type = rough
     roughness = {{ epsilon }}
     gravity = '0 0 -9.8'
     block = r430_2150
     outputs = exodus
   [../]
   [./lat_r430_2150]
     type = MoskitoLatHeat_Inc_Formation_1p
     temperature_inner = T
     outer_diameters = '{{ diameter }} 0.1778 0.2168 0.2450 0.3179'
     conductivities = '43.75 0.7 43.75 0.7'
     formation_density = {{ rock_density }}
     formation_heat_capacity = {{ rock_cp }}
     formation_temperature_function = grad_func
     convective_thermal_resistance = true
     nondimensional_time_function = Hasan_Kabir_2012
     outputs = exodus
     block = r430_2150
     formation_thermal_conductivity = {{ rock_k }}
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
      well_diameter = {{ diameter }}
      roughness_type = rough
      roughness = {{ epsilon }}
      gravity = '0 0 -9.8'
      block = r0_430
      outputs = exodus
    [../]
    [./lat_r0_430]
      type = MoskitoLatHeat_Inc_Formation_1p
      temperature_inner = T
      outer_diameters = '{{ diameter }} 0.1778 0.2168 0.2450 0.3179 0.3397 0.4'
      conductivities = '43.75 0.3 43.75 0.7 43.75 0.7'
      formation_density = {{ rock_density }}
      formation_heat_capacity = {{ rock_cp }}
      formation_temperature_function = grad_func
      convective_thermal_resistance = true
      nondimensional_time_function = Hasan_Kabir_2012
      outputs = exodus
      block = r0_430
      formation_thermal_conductivity = {{ rock_k }}
    [../]
[]

[Functions]
  [./dts]
    type = PiecewiseLinear
    x = '0 20000  315360  31536000 630720000'
    y = '2000.0 2000.0 86400 3153600 1000000'
  [../]
  [./conc]
    type = ParsedFunction
    value = 'if(t>16000 & t<73000, 1.72, 0.25)'
  [../]
  [./grad_func]
    type = ParsedFunction
    value ='{{ surface_temperature }} - z * {{ thermal_gradient }}'
  [../]
  [./conductivity_gradient]
    type = ParsedFunction
    value ='1.5 - z * 0.00075'
  [../]
  [./csv_data]
    type = MoskitoConstant
    data_file = pc.csv
    direction = RIGHT
    format = columns
    execute_on = TIMESTEP_BEGIN
  [../]
  [./csv_datat]
    type = MoskitoConstant
    data_file = tc.csv
    direction = RIGHT_INCLUSIVE
    format = columns
    execute_on = TIMESTEP_BEGIN
  [../]
[]

[BCs]
#  [./pleft]
#    type = DirichletBC
#    variable = p
#    boundary = inlet
#    value = 10000000
#  [../]
  [./pleft]
    type = FunctionDirichletBC
    variable = p
    boundary = inlet
    function = csv_data
  [../]
  [./qleft]
    type = FunctionDirichletBC
    variable = q
    boundary = inlet
    function = '0.01'
  [../]
#  [./tleft]
#    type = DirichletBC
#    variable = T
#    boundary = inlet
#    value = 300.0
#  [../]
   [./tleft]
    type = FunctionDirichletBC
    variable = T
    boundary = inlet
    function = csv_datat
  [../]
[]

[Variables]
  [./T]
    [./InitialCondition]
      type = FunctionIC
      function = '{{ surface_temperature }} - z * {{ thermal_gradient }}'
      variable = T
    [../]
  [../]
  [./p]
    [./InitialCondition]
      type = FunctionIC
      function = '8000000 -10000 * z'
      variable = p
    [../]
  [../]
  [./q]
    initial_condition = 0.01
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
  end_time = {{ end_time }}
  l_max_its = 100
  solve_type = NEWTON
  steady_state_detection = false
  [./TimeStepper]
    type = FunctionDT
    function = dts
  [../]
[]

[Outputs]
  file_base = new_out
  exodus = true
  csv = true
[]

[Postprocessors]
  [belphip]
    type = PointValue
    variable = p
    point = '5000 0 0'
    outputs = 'csv'
  []
  [belphit]
    type = PointValue
    variable = T
    point = '5000 0 0'
    outputs = 'csv'
  []
  [belphiq]
    type = PointValue
    variable = q
    point = '5000 0 0'
    outputs = 'csv'
  []
[]
