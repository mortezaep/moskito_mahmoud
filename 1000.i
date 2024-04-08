# Validation of Heat exchange using the paper from Satman & Tureyen 2016
# Actually ramey analytical solution
[Mesh]
  type = FileMesh
   file = swm.msh
[]

[UserObjects]
  [./viscosity]
    type = MoskitoViscosityWaterVogel
  [../]
  [./eos]
    type = MoskitoEOS1P_Brine
  [../]
  [./annulus]
    type = MoskitoAnnulus
    density = 1000.0
    viscosity = 0.001
    thermal_expansion = 2.0e-4
    thermal_conductivity = 0.6
    heat_capacity = 4000.0
    emissivity_inner = 0.9
    emissivity_outer = 0.9
    convective_heat_method = Churchill
  [../]
[]

[Functions]
  [./dts]
    type = PiecewiseLinear
    #x = '0 7200  72000.0 943000  315360000 3153600000'
    #y = '100 1000.0 1000.0 360000 3600000 3600000'
     x = '0     1576800000'
     y = '3600  8640000'
  [../]
  [./grad_func]
    type = ParsedFunction
    value ='281.51 - z * 0.03598'
  [../]
  [./conductivity_gradient]
    type = ParsedFunction
    value = 'if(z>-2000, 3.1, 3.1)'
  [../]
[]

[Materials]
  [./area1]
    type = MoskitoFluidWell_1p1c
    temperature = T
    pressure = p
    flowrate = q
    well_direction = -z
    well_diameter = 0.2159
    eos_uo = eos
    viscosity_uo = viscosity
    roughness_type = smooth
    gravity = '0 0 -9.8'
    well_type = -1
    block = '3491_3926'
  [../]
  [./area2]
    type = MoskitoFluidWell_1p1c
    temperature = T
    pressure = p
    flowrate = q
    well_direction = -z
    well_diameter = 0.244475
    eos_uo = eos
    viscosity_uo = viscosity
    roughness_type = smooth
    gravity = '0 0 -9.8'
    well_type = -1
    block = '2179_3491'
  [../]
  [./area3]
    type = MoskitoFluidWell_1p1c
    temperature = T
    pressure = p
    flowrate = q
    well_direction = -z
    well_diameter = 0.244475
    eos_uo = eos
    viscosity_uo = viscosity
    roughness_type = smooth
    gravity = '0 0 -9.8'
    well_type = -1
    block = '2072_2179'
  [../]
  [./area4]
    type = MoskitoFluidWell_1p1c
    temperature = T
    pressure = p
    flowrate = q
    well_direction = -z
    well_diameter = 0.339725
    eos_uo = eos
    viscosity_uo = viscosity
    roughness_type = smooth
    gravity = '0 0 -9.8'
    well_type = -1
    block = '953_2072'
  [../]
  [./area5]
    type = MoskitoFluidWell_1p1c
    temperature = T
    pressure = p
    flowrate = q
    well_direction = -z
    well_diameter = 0.339725
    eos_uo = eos
    viscosity_uo = viscosity
    roughness_type = smooth
    gravity = '0 0 -9.8'
    well_type = -1
    block = '892_953'
  [../]
  [./area6]
    type = MoskitoFluidWell_1p1c
    temperature = T
    pressure = p
    flowrate = q
    well_direction = -z
    well_diameter = 0.244475
    eos_uo = eos
    viscosity_uo = viscosity
    roughness_type = smooth
    gravity = '0 0 -9.8'
    well_type = -1
    block = '802_892'
  [../]
  [./area7]
    type = MoskitoFluidWell_1p1c
    temperature = T
    pressure = p
    flowrate = q
    well_direction = -z
    well_diameter = 0.244475
    eos_uo = eos
    viscosity_uo = viscosity
    roughness_type = smooth
    gravity = '0 0 -9.8'
    well_type = -1
    block = '48_802'
  [../]
  [./area8]
    type = MoskitoFluidWell_1p1c
    temperature = T
    pressure = p
    flowrate = q
    well_direction = -z
    well_diameter = 0.244475
    eos_uo = eos
    viscosity_uo = viscosity
    roughness_type = smooth
    gravity = '0 0 -9.8'
    well_type = -1
    block = '0_48'
  [../]
    [./Lateral1]
      type = MoskitoLatHeat_Inc_Formation_1p
       temperature_inner = T
       outer_diameters = '0.2159 0.21591'
       conductivities = '3.1'
       convective_thermal_resistance = true
       # Rock parameters
       formation_density = 2500
       formation_thermal_conductivity = 3.1
       formation_heat_capacity = 900
       # Configuration of material
       formation_temperature_function = grad_func
       nondimensional_time_function = Ramey_1981_BF
       output_properties = 'total_thermal_resistivity'
       outputs = exodus
       block = '3491_3926'
     [../]
     [./Lateral2]
       type = MoskitoLatHeat_Inc_Formation_1p
        temperature_inner = T
        outer_diameters = '0.244475 0.254475 0.31115'
        conductivities = '45.7 1.1'
        convective_thermal_resistance = true
        # Rock parameters
        formation_density = 2500
        formation_thermal_conductivity = 3.1
        formation_heat_capacity = 900
        # Configuration of material
        formation_temperature_function = grad_func
        nondimensional_time_function = Ramey_1981_BF
        output_properties = 'total_thermal_resistivity'
        outputs = exodus
        block = '2179_3491'
      [../]
      [./Lateral3]
        type = MoskitoLatHeat_Inc_Formation_1p
         temperature_inner = T
         outer_diameters = '0.244475 0.254475 0.339725 0.350625 0.4064'
         conductivities = '45.7 1.1 45.7 1.1'
         convective_thermal_resistance = true
         # Rock parameters
         formation_density = 2500
         formation_thermal_conductivity = 3.1
         formation_heat_capacity = 900
         # Configuration of material
         formation_temperature_function = grad_func
         nondimensional_time_function = Ramey_1981_BF
         output_properties = 'total_thermal_resistivity'
         outputs = exodus
         block = '2072_2179'
       [../]
       [./Lateral4]
         type = MoskitoLatHeat_Inc_Formation_1p
          temperature_inner = T
          outer_diameters = '0.339725 0.350625 0.4064'
          conductivities = '45.7 1.1'
          convective_thermal_resistance = true
          # Rock parameters
          formation_density = 2500
          formation_thermal_conductivity = 3.1
          formation_heat_capacity = 900
          # Configuration of material
          formation_temperature_function = grad_func
          nondimensional_time_function = Ramey_1981_BF
          output_properties = 'total_thermal_resistivity'
          outputs = exodus
          block = '953_2072'
        [../]
        [./Lateral5]
          type = MoskitoLatHeat_Inc_Formation_1p
           temperature_inner = T
           outer_diameters = '0.339725 0.350625 0.473075 0.483075 0.5842'
           conductivities = '45.7 1.1 45.7 0.02'
           convective_thermal_resistance = true
           # Rock parameters
           formation_density = 2500
           formation_thermal_conductivity = 3.1
           formation_heat_capacity = 900
           # Configuration of material
           formation_temperature_function = grad_func
           nondimensional_time_function = Ramey_1981_BF
           output_properties = 'total_thermal_resistivity'
           outputs = exodus
           block = '892_953'
         [../]
         [./Lateral6]
           type = MoskitoLatHeat_Inc_Formation_1p
            temperature_inner = T
            outer_diameters = '0.244475 0.254475 0.339725 0.351925 0.473075 0.483075 0.5842'
            conductivities = '45.7 0.0 45.7 1.1 45.7 0.02'
            convective_thermal_resistance = true
            # Rock parameters
            formation_density = 2500
            formation_thermal_conductivity = 3.1
            formation_heat_capacity = 900
            # Configuration of material
            formation_temperature_function = grad_func
            nondimensional_time_function = Ramey_1981_BF
            output_properties = 'total_thermal_resistivity'
            outputs = exodus
            block = '802_892'
            annulus_uo = annulus
          [../]
          [./Lateral7]
            type = MoskitoLatHeat_Inc_Formation_1p
             temperature_inner = T
             outer_diameters = '0.244475 0.254475 0.339725 0.351925 0.473075 0.483075 0.5842'
             conductivities = '45.7 0.6 45.7 0.0 45.7 0.02'
             convective_thermal_resistance = true
             # Rock parameters
             formation_density = 2500
             formation_thermal_conductivity = 3.1
             formation_heat_capacity = 900
             # Configuration of material
             formation_temperature_function = grad_func
             nondimensional_time_function = Ramey_1981_BF
             output_properties = 'total_thermal_resistivity'
             outputs = exodus
             block = '48_802'
             annulus_uo = annulus
           [../]
           [./Lateral8]
             type = MoskitoLatHeat_Inc_Formation_1p
              temperature_inner = T
              outer_diameters = '0.244475 0.254475 0.339725 0.351925 0.473075 0.483075 0.6096 0.6196 0.6604'
              conductivities = '45.7 0.6 45.7 0.0 45.7 0.02 45.7 1.1'
              convective_thermal_resistance = true
              # Rock parameters
              formation_density = 2500
              formation_thermal_conductivity = 3.1
              formation_heat_capacity = 900
              # Configuration of material
              formation_temperature_function = grad_func
              nondimensional_time_function = Ramey_1981_BF
              output_properties = 'total_thermal_resistivity'
              outputs = exodus
              block = '0_48'
              annulus_uo = annulus
            [../]
[]

[Variables]
  [./T]
    [./InitialCondition]
      type = FunctionIC
      function = '283.15 - z * 0.03'
      #function = '281.51 - z * 0.03598'
      variable = T
    [../]
  [../]
  [./p]
    [./InitialCondition]
      type = FunctionIC
      function = '-10000 * z + 1500000'
      variable = p
    [../]
  [../]
  [./q]
    initial_condition = 0.115
  [../]
[]

[BCs]
  [./Tbc]
    type = DirichletBC
    variable = T
    boundary = bot
    value =  414.35
  [../]
  [./pleft]
    type = DirichletBC
    variable = p
    boundary = top
    value = 1500000
  [../]
  [./qleft]
    type = DirichletBC
    variable = q
    boundary = bot
    value = 0.115
  [../]
[]

[Kernels]
  [./Tkernel]
    type = MoskitoEnergy_1p1c
    variable = T
    pressure = p
    flowrate = q
  [../]
  [./Ttkernel]
    type = MoskitoTimeEnergy_1p1c
    variable = T
    pressure = p
    flowrate = q
  [../]
  [./heat]
    type = MoskitoLatHeatIncFormation_1p
    variable = T
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
[]

[Preconditioning]
  active = p3
  [./p3]
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-pc_type -ksp_type -sub_pc_type -snes_type -snes_linesearch_type -pc_asm_overlap -sub_pc_factor_shift_type -ksp_gmres_restart'
    petsc_options_value = 'asm gmres lu newtonls basic 2 NONZERO 51'
  [../]
[]

[Executioner]
  type = Transient
  start_time = 0
  end_time = 1576800000
  l_max_its = 100
  solve_type = NEWTON
  nl_abs_tol = 1e-5
  # nl_abs_tol = 1e5
  [./TimeStepper]
    type = FunctionDT
    function = dts
  [../]
[]

[Outputs]
  exodus = true
  print_linear_residuals = false
[]
