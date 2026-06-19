import numpy as np

# use log scaling when log difference between values at the lower bound of the range (e.g., 1 and 10) is as physically significant as at the upper bond (e.g., 10000 and 100,000).
# use linear scaling when unit difference between values at the lower bound of the range (e.g., 0.1 and 0.2) is as physically significant as at the upper bond (e.g., 15.9 and 16.0). Can't use zeroes in range.

# target = base * ratio, e.g. L = D * L/D or Dp = D * Dp/D

input_params = [
# packed bed design
  {'name': 'D',            'logScaling': False,  'type': 'independent',                          'range': (0.1, 5),             'units': 'm',        'desc': 'Packed bed equivalent diameter'},
  {'name': 'L',            'logScaling': True,   'type': 'dependent',    'base': 'D',            'ratio_range': (1, 500),       'units': 'm',        'desc': 'Packed bed length'},
  {'name': 'D_p',          'logScaling': True,   'type': 'dependent',    'base': 'D',            'ratio_range': (1/500, 1/1.5), 'units': 'm',        'desc': 'Particle equivalent diameter'},
  {'name': 'phi_s',        'logScaling': False,  'type': 'independent',                          'range': (0.5, 1),             'units': None,       'desc': 'Particle sphericity'},
  {'name': 'ε',            'logScaling': False,  'type': 'independent',                          'range': (0.36, 0.9),          'units': None,       'desc': 'Dry packed bed void fraction'},

# gaseous phase (air)
  {'name': 'U_s_G_in',      'logScaling': True,   'type': 'independent',                          'range': (0.01, 5),           'units': 'm/s',      'desc': 'Superficial velocity of the gaseous phase at the inlet'},
]

output_params = [
  {'name': 'ΔP',            'units': 'Pa',       'desc': 'Pressure drop'},
]

scaling_lookup = {p['name']: p['logScaling'] for p in input_params}
scaling_lookup.update({
  'Re_m': True,   # Reynolds number spans orders of magnitude -> Log
  'h_L': False,   # Liquid holdup is bounded between 0 and 1 -> Linear
})

############################################################################################
# Source Models
#
# Model Types:
# - analytical     — pure math without any ackward coeffs
# - semi-empirical — elegant generalized correlations but with ackward coeffs
# - empirical      — meaningless correlations with tight validity ranges
############################################################################################

src_models = [
  # 'inContextLinear' is a list of input params with unit difference between values at the lower bound of the range (e.g., 0.1 and 0.2) is as physically significant as at the upper bond (e.g., 15.9 and 16.0).
  # 'inContextLog' is a list of input params with log difference between values at the lower bound of the range (e.g., 1 and 10) is as physically significant as at the upper bond (e.g., 10000 and 100,000).
  # 'outContext' is a list of output params.
  # Re_p context includes U_s_G_in, D_p, ρ_G_avg, and U_s_G_in.

  ##### PRESSURE DROP ##########################################################################################################################
  {'name':'test',             'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'carman1937',       'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'D', 'phi_s', 'ε'],        'outContext':['ΔP']},
  {'name':'morcom1946',       'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'rose1949',         'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'leva1949',         'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'ergun1952',        'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'wentz1963',        'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'kurten1966',       'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'handley1968',      'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'mehta1969',        'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'D', 'ε'],                 'outContext':['ΔP']},
  {'name':'hicks1970',        'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'tallmadge1970',    'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'reichelt1972',     'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'D', 'ε'],                 'outContext':['ΔP']},
  {'name':'kuo1978',          'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'macdonald1979',    'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'kta1981',          'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'jones1983',        'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'foscolo1983',      'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'fahien1983',       'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'meyer1985',        'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'phi_s', 'ε'],             'outContext':['ΔP']},
  {'name':'paterson1986',     'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'D', 'ε'],                 'outContext':['ΔP']},
  {'name':'stichlmair1989',   'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε', 'h_L'],               'outContext':['ΔP']},
  {'name':'watanabe1989',     'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'comiti1989',       'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'D', 'ε'],                 'outContext':['ΔP']},
  {'name':'fand1990',         'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'D', 'ε'],                 'outContext':['ΔP']},
  {'name':'foumeny1993',      'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'D', 'ε'],                 'outContext':['ΔP']},
  {'name':'lee1994',          'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'liu1994',          'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'D', 'ε'],                 'outContext':['ΔP']},
  {'name':'hayes1995',        'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'avontuur1996',     'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'critoph1996',      'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'oneill1997',       'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'D', 'ε'],                 'outContext':['ΔP']},
  {'name':'raichura1999',     'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'D', 'ε'],                 'outContext':['ΔP']},
  {'name':'eisfeld2001',      'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'D', 'phi_s', 'ε'],        'outContext':['ΔP']},
  {'name':'yu2002',           'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'felice2004',       'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'D', 'ε'],                 'outContext':['ΔP']},
  {'name':'nemec2005',        'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'phi_s', 'ε'],             'outContext':['ΔP']},
  {'name':'montillet2007',    'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'D', 'ε'],                 'outContext':['ΔP']},
  {'name':'carpinlioglu2008', 'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'ozahi2008',        'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'phi_s', 'ε'],             'outContext':['ΔP']},
  {'name':'cheng2011',        'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'D', 'ε'],                 'outContext':['ΔP']},
  {'name':'harrison2013',     'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'D', 'ε'],                 'outContext':['ΔP']},
  {'name':'erdim2015',        'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'guo2017',          'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'D', 'ε'],                 'outContext':['ΔP']},
  {'name':'seckendorff2020',  'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'cheng2021',        'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'D', 'ε'],                 'outContext':['ΔP']},
  {'name':'reger2023',        'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'dixon2023',        'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'D', 'ε'],                 'outContext':['ΔP']},
  {'name':'wu2025',           'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'ε'],                      'outContext':['ΔP']},
  {'name':'xu2026',           'src_type':'semi-empirical', 'inContext':['Re_p', 'L', 'phi_s', 'ε'],             'outContext':['ΔP']},

  ##### HEAT & MASS TRANSFER ###################################################################################################################

]

a_rules = [
  # 'y' is an initial output param (e.g., ΔP)
  # 'input' is a function to define input param (e.g., U²ρ/D_p)
  # 'output' is a function to transform initial output param (e.g., ΔP → ΔP/L)
  # 'expected_corr' is a sign of correlation gradient ∂Y/∂X (`pos` for positive correlation and `neg` for negative correlation)
  # 'control_param' is an input/output param used to define the validity window of the rule
  # 'control_range' is a validity range of 'control_param'
  # 'context' is a list of input/output params covered by the rule. Re_p context includes U_s_G_in, D_p, ρ_G_avg, and U_s_G_in.

  {
    'name': 'Laminar Stress Test', ## turb and lam regimes are not merged into a single rule because neighbourhood definition doesn't respect jumps created by np.where switch.
    'y': 'ΔP',
    'input': lambda dm, **kwargs: dm['U_s_G_in'] * dm['μ_G_avg'] / dm['D_p']**2,
    'output': lambda dm, y: y / dm['L'],
    'expected_corr': 'pos',
    'control_param': 'Re_m',
    'control_range': (None, 1000),
    'context':['Re_p', 'L', 'ΔP']
  },
  {
    'name': 'Turbulent Stress Test',
    'y': 'ΔP',
    'input': lambda dm, **kwargs: dm['U_s_G_in']**2 * dm['ρ_G_avg'] / dm['D_p'],
    'output': lambda dm, y: y / dm['L'],
    'expected_corr': 'pos',
    'control_param': 'Re_m',
    'control_range': (10, None),
    'context':['Re_p', 'L', 'ΔP']
  },
]

b_rules = [
  
]