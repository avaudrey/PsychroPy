#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Python package dedicated to calculations dealing with moist air and HVAC
application, based on the CoolProp package (http://www.coolprop.org/).

Author: alexandre.vaudrey@gmail.com
"""

#----------------------------------------------------------------------------
#   Copyright (C) 2018 <Alexandre Vaudrey>                                  |
#                                                                           |
#   This program is free software: you can redistribute it and/or modify    |
#   it under the terms of the GNU General Public License as published by    |
#   the Free Software Foundation, either version 3 of the License, or       |
#   (at your option) any later version.                                     |
#                                                                           |
#   This program is distributed in the hope that it will be useful,         |
#   but WITHOUT ANY WARRANTY; without even the implied warranty of          |
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           |
#   GNU General Public License for more details.                            |
#                                                                           |
#   You should have received a copy of the GNU General Public License       |
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.   |
#---------------------------------------------------------------------------|

# Properties of humid air
from CoolProp.CoolProp import HAPropsSI
# Value of one atmosphere
ATM = 1.01325

#TODO : find how to use alas for main properties of air
# specific humidity = humidity ratio = W
# relative humidity = R
# specific enthalpy = H
# specific entropy = S
# Wet bulb temperature = B
# Dew point temperature = D

class MoistAir(object):
    """
    Object containing the state of a given moist air. Temperature is expressed
    in °C and pressure in bar. Although three parameters, as e.g. temperature
    (T), relative humidity (R) and pressure (P) are required to identify the
    specific state of any moist air, pressure will be implicitly considered as
    equal to its default value, i.e. 1 atm, until a different value is
    specifically entered as an attribute. A couple of parameter as e.g. RT, RW,
    TW, can then be used to enter the specific state of the air (as it is used
    for example in CANTERA) with:
        R = relative humidity
        T = dry bulb temperature
        W = specific humidity or humidity ratio
        B = wet bulb temperature
        D = dew point temperature
        H = specific enthalpy
    See the README file for further information.
    """
    def __init__(self):
        # Name of the state
        self.name = 'State 1'
        # Pressure
        self.pressure = 1*ATM
        # Temperature
        self.temperature = 20.
        # Relative humidity
        self.relative_humidity = 0.5
        # Specific humidity
        self.specific_humidity = HAPropsSI('W',
                                           'P', self.pressure*1e5,
                                           'T', self.temperature+273.15,
                                           'R', self.relative_humidity)
    # TODO : create a function that can read the name of any attribute entered
    # as PWT, PRT, BTR,... and so on, identify each parameter from a dictionary,
    # and attribute the joined value to the corresponding parameter.

if __name__ == '__main__':
    inside = MoistAir()
