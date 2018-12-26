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
import sys
# Value of one atmosphere
ATM = 1.01325

#TODO : find how to use alas for main properties of air
# specific humidity = humidity ratio = W
# relative humidity = R
# specific enthalpy = H
# specific entropy = S
# Wet bulb temperature = B
# Dew point temperature = D

# Dictionnary of properties
properties = {'R':'relative_humidity',
              'W':'specific_humidity',
              'T':'temperature',
              'B':'wet_bulb_temperature',
              'D':'dew_point_temperature',
              'H':'specific_enthalpy',
              'P':'pressure'}

single_parameter_set_message = "At least two complementary properties of moist air are required to unambiguously identify the latter, any sole parameter can then not be used for so."

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
        self._specific_humidity = HAPropsSI('W',
                                            'P', self.pressure*1e5,
                                            'T', self.temperature+273.15,
                                            'R', self.relative_humidity)
    # At least two values being required to identify unambiguously the state of
    # a given moist air, any new moist air must be entered thanks to a couple of
    # properties, as e.g. the dry and the wet bulb temperature or the dry
    # temperature and the relative humidity. Both single property of moist air
    # can then be get but not set directly.
    @property
    def specific_humidity(self):
        """
        Specific humidity, aka humidity ratio, defined as the ratio of the mass
        of water vapour on the dry air one it is mixed with. Dimensionless.
        """
        return self._specific_humidity
    @specific_humidity.setter
    def specific_humidity(self, value):
        print(single_parameter_set_message)
        pass
    # Dry temperature and relative humidity as entered parameters
    @property
    def TR(self):
        """
        Dry temperature and relative humidity as entered properties of the moist
        air.
        """
        return self.temperature, self.relative_humidity
    @TR.setter
    def TR(self, values):
        """
        Dry temperature and relative humidity as known properties of moist air.
        """
        # Set of the temperature value
        self.temperature = values[0]
        # Check if the entered value of relative humidity is between 0 and 1
        if (values[1] < 0):
            self.relative_humidity = 0.0
            print('Relative humidity is a number between 0% and 100%!')
        elif (values[1] > 1):
            self.relative_humidity = 1.0
            print('Relative humidity is a number between 0% and 100%!')
        else:
            self.relative_humidity = values[1]
        # Calculation of the corresponding specific humidity value
        self._specific_humidity = HAPropsSI('W',
                                            'P', self.pressure*1e5,
                                            'T', self.temperature+273.15,
                                            'R', self.relative_humidity)
        pass
    @property
    def RT(self):
        """
        Relative humidity and dry temperature as entered properties of the moist
        air.
        """
        return self.relative_humidity, self.temperature
    @RT.setter
    def RT(self, values):
        """
        Relative humidity and dry temperature as entered properties of the moist
        air.
        """
        self.TR = values[::-1]
        pass

if __name__ == '__main__':
    inside = MoistAir()
    inside.TR = 0.0, 0.8
    print(inside.specific_humidity)
