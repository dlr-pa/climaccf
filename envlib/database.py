#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# temperature threshold for foramtion of contrails
cont_temp_thr = 235

# relative hunmidity threshold for presistancy of formed contrails
cont_rh_thr = 1

# threshold for determining climate hotspots
threshold = 1e-14

# Forcing efficacies (lee et al. 2021)
efficacy = {'CH4': 1.18, 'O3': 1.37, 'H2O': 1, 'Cont.': 0.42, 'CO2': 1}

# Typical transatlantic values (TTV)
emission_index = {'Cont.': 0.16, 'NOx': 13 * 1e-3}

# Pulse to future
P2F = {'CH4': 9.61, 'O3': 12.66, 'H2O': 12.66, 'Cont.': 12.3, 'CO2': 9.28}
