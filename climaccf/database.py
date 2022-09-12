#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Fixed threshold for determining climate hotspots
threshold = 1e-14

# Forcing efficacies after lee et al. 2021
efficacy = {'CH4': 1.18, 'O3': 1.37, 'H2O': 1, 'Cont.': 0.42, 'CO2': 1}

# Typical transatlantic fleet mean values (TTV) of NOx emission index and flown distance
emission_index = {'Cont.': 0.16, 'NOx': 13 * 1e-3}

# Metric conversion factors form pulse emission to future emission scenario over different time horizons (see Dietmueller et al. 22)
P20_F20 = {'CH4': 10.8, 'O3': 14.5, 'H2O': 14.5, 'Cont.': 13.6, 'CO2': 9.4}
P20_F50 = {'CH4': 42.5, 'O3': 34.1, 'H2O': 34.1, 'Cont.': 30.1, 'CO2': 44.0}
P20_F100 = {'CH4': 98.2, 'O3': 58.3, 'H2O': 58.3, 'Cont.': 48.9, 'CO2': 125.0}


# Scaling factors for the first complete and consistent set of aCCFs  (aCCF-V1.0)
aCCF_V1_0 = {'CH4':2.03 , 'O3':1.97 , 'H2O':1.92 , 'Cont.': 1, 'CO2': 1}

# Educated guess factors after Matthes et al. 2022 (aCCF-V1.1)
aCCF_V1_1 = {'CH4': 35, 'O3': 11, 'H2O': 3, 'Cont.': 3, 'CO2': 1} 

