[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = -1
  xmax = 1
  ymin = -1
  ymax = 1
  nx = 10
  ny = 10
[]

[Variables]
  [./u]
    [./InitialCondition]
      type = ConstantIC
      value = 0
    [../]
  [../]
[]

[Functions]
  # [./csv_reader]
  #   type = PiecewiseBilinear
  #   data_file = bcs.csv
  #   yaxis = 2
  # [../]
  [./csv_data]
    type = MoskitoConstant
    data_file = pcres.csv
    direction = RIGHT
    format = columns
    execute_on = TIMESTEP_BEGIN
  [../]
[]

[Kernels]
  [./ie]
    type = TimeDerivative
    variable = u
  [../]
  [./diff]
    type = Diffusion
    variable = u
  [../]
[]

[BCs]
  [./all]
    type = FunctionDirichletBC
    variable = u
    boundary = left
    function = csv_data
  [../]
  [./right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  [../]
[]

[Postprocessors]
  [belphires]
    type = PointValue
    variable = u
    point = '-1 0 0'
  []
[]

[Executioner]
  type = Transient
  solve_type = Newton
  start_time = 0.0
  num_steps = 20
  dt = 0.1
[]

[Outputs]
  execute_on = 'timestep_end'
  file_base = outres_transient
  exodus = true
  csv = true
[]

