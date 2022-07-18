#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# temperature threshold for foramtion of persistent contrails [in K]
cont_temp_thr = 235

# relative hunmidity threshold for formation of persistent contrails
cont_rh_thr = 1

# fixed threshold for determining climate hotspots
threshold = 1e-14

# Forcing efficacies after lee et al. 2021
efficacy = {'CH4': 1.18, 'O3': 1.37, 'H2O': 1, 'Cont.': 0.42, 'CO2': 1}

# Typical transatlantic fleet mean values (TTV) of NOx emission index and flown distance
emission_index = {'Cont.': 0.16, 'NOx': 13 * 1e-3}

# Metric conversion factors form pulse emission to future emission scenario over differnt time horizons (see Dietmueller et al. 22)

P20_F20 = {'CH4': 9.83, 'O3': 14.47, 'H2O': 14.47, 'Cont.': 13.42, 'CO2': 9.28}
P20_F50 = {'CH4': 41.48, 'O3': 34.13, 'H2O': 34.13, 'Cont.': 29.46, 'CO2': 43.92}
P20_F100 = {'CH4': 97.37, 'O3': 58.34, 'H2O': 58.34, 'Cont.': 47.31, 'CO2': 124.93}

# educated guess factors after Matthes et al. 2022 (aCCF-V1.1)
educated_guess_v1 = {'CH4': 35, 'O3': 11, 'H2O': 3, 'Cont.': 3, 'CO2': 1}
