#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Python package dedicated to calculations dealing with moist air and HVAC
application, mostly based on the CoolProp package (http://www.coolprop.org/).

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

import numpy as np
# Properties of humid air
from CoolProp.CoolProp import HAPropsSI
# Value of one atmosphere in bar
ATM = 1.01325
# Specific heats at constant pressure of dry air and water vapor, in [J/(kg.K)]
DRY_AIR_CP, WATER_VAPOR_CP = 1004., 1805.
# Specific gas constants of dry air and water vapor, in [J/(kg.K)]
DRY_AIR_R, WATER_VAPOR_R = 287., 462.
# Ratio of the water molar mass on the dry air one
ALPHAW = DRY_AIR_R/WATER_VAPOR_R
# Water specific enthalpy of vaporization, in [J/(kg.K)]
WATER_LW = 2501e+3

single_parameter_set_message = "At least two complementary properties of moist air are required to unambiguously identify the latter, any sole parameter can then not be used for so."

class MoistAir(object):
    """
    Object containing the state of a given moist air. Temperature is expressed
    in °C, pressure is expressed in bar and specific enthalpy in kJ/kg. Although
    three parameters, as e.g. temperature (T), relative humidity (R) and
    pressure (P) are theoretically required to identify the specific state of
    any moist air, pressure will be implicitly considered as equal to its
    default value, i.e. 1 atm, until a different value is specifically entered
    as attribute. A couple of parameter as e.g. RT, RW, TW, must then be used to
    enter the specific state of the air (as it is used for example in CANTERA)
    with:
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
        # When the class is involved in computations dealing with a large amount
        # of data, the below attribute can be switched to True in order to
        # accelerate the calculation is avoiding as far as possible the use of
        # the CoolProp package. Obtained results are usually a little bit less
        # accurate than the ones obtained by Coolprop
        self.fast_computation = False
        # Pressure, in bar
        self.pressure = 1*ATM
        # Temperature, in °C
        self._temperature = 20.
        # Relative humidity, dimensionless
        self._relative_humidity = 0.5
        # Specific humidity, dimensionless
        self._specific_humidity = HAPropsSI('W',
                                            'P', self.pressure*1e5,
                                            'T', self.temperature+273.15,
                                            'R', self.relative_humidity)
        # Wet bulb temperature, in °C
        self._wet_bulb_temperature = HAPropsSI('B',
                                               'P', self.pressure*1e5,
                                               'T', self.temperature+273.15,
                                               'R', self.relative_humidity)\
                -273.15
        # Dew point temperature, in °C
        self._dew_point_temperature = HAPropsSI('D',
                                                'P', self.pressure*1e5,
                                                'T', self.temperature+273.15,
                                                'R', self.relative_humidity)\
                -273.15
        # Specific enthalpy, in kJ/kg
        self._specific_enthalpy = HAPropsSI('H',
                                            'P', self.pressure*1e5,
                                            'T', self.temperature+273.15,
                                            'R', self.relative_humidity)*1e-3
        # Specific volume, in m3/kg 
        self._specific_volume = HAPropsSI('V',
                                          'P', self.pressure*1e5,
                                          'T', self.temperature+273.15,
                                          'R', self.relative_humidity)
    # Static functions ========================================================
    @staticmethod
    def __equilibrium_vapor_pressure(temperature):
        """
        Saturation or equilibrium pressure of water, in Pa.
        """
        # TODO: Find the original reference where this equation actually comes
        # from.
        # [ref?] equation that gives a value of the equilibrium pressure of
        # water with a maximum difference of 0.6% (when compared with CoolProp)
        # on a temperature range such as 5°C < T < 100°C. This equation is used
        # when fast computation are required.
        return 1e5*np.exp(11.78*(temperature-99.64)/(temperature+230))
    @staticmethod
    def __logmeantemperature(temp1, temp2):
        """
        Log mean temperature of temperature temp1 and temp2, in K, with both
        temperature entered in °C.
        """
        # This function is used for example in the calculation of the specific
        # thermal exergy of moist air.
        if temp1 == temp2:
            lmtemp = temp1+273.15
        else:
            lmtemp = (temp1-temp2)/np.log((temp1+273.15)/(temp2+273.15))
        return lmtemp
    # ==== Single properties of moist air =====================================
    # At least two values being required to identify unambiguously the state of
    # a given moist air, any new moist air must be entered thanks to a couple of
    # properties, as e.g. the dry and the wet bulb temperature or the dry
    # temperature and the relative humidity. Both single property of moist air
    # can then be get but not set directly.
    # ---- Temperature --------------------------------------------------------
    @property
    def temperature(self):
        """
        Usual (dry bulb) temperature, in °C.
        """
        return self._temperature
    @temperature.setter
    def temperature(self, value):
        print(single_parameter_set_message)
    @temperature.deleter
    def temperature(self):
        raise AttributeError("Can't delete attribute")
    # ---- Relative humidity --------------------------------------------------
    @property
    def relative_humidity(self):
        """
        Relative humidity, dimensionless, defined as the ratio of the actual
        water vapour partial pressure on the maximum one equal to the
        equilibrium pressure.
        """
        return self._relative_humidity
    @relative_humidity.setter
    def relative_humidity(self, value):
        print(single_parameter_set_message)
    @relative_humidity.deleter
    def relative_humidity(self):
        raise AttributeError("Can't delete attribute")
    # ---- Specific humidity --------------------------------------------------
    @property
    def specific_humidity(self):
        """
        Specific humidity, aka humidity ratio, defined as the ratio of the mass
        of water vapour on the dry air one it is mixed with. Dimensionless.
        """
        return self._specific_humidity
    @specific_humidity.setter
    def specific_humidity(self):
        print(single_parameter_set_message)
    @specific_humidity.deleter
    def specific_humidity(self):
        raise AttributeError("Can't delete attribute")
    # ---- Wet bulb temperature -----------------------------------------------
    @property
    def wet_bulb_temperature(self):
        """
        Wet bulb temperature in °C, se the one of moist air if cooled down to
        saturation through an adiabatic process, i.e. with a constant specific
        enthalpy.
        """
        return self._wet_bulb_temperature
    @wet_bulb_temperature.setter
    def wet_bulb_temperature(self, value):
        print(single_parameter_set_message)
    @wet_bulb_temperature.deleter
    def wet_bulb_temperature(self, value):
        raise AttributeError("Can't delete attribute")
    # ---- Dew point temperature ----------------------------------------------
    @property
    def dew_point_temperature(self):
        """
        Dew point temperature in °C, so the one of moist air if cooled down to
        saturation with a constant specific humidity.
        """
        return self._dew_point_temperature
    @dew_point_temperature.setter
    def dew_point_temperature(self, value):
        print(single_parameter_set_message)
    @dew_point_temperature.deleter
    def dew_point_temperature(self, value):
        raise AttributeError("Can't delete attribute")
    # ---- Specific enthalpy --------------------------------------------------
    @property
    def specific_enthalpy(self):
        """
        Specific enthalpy of moist air, in kJ/kg.
        """
        return self._specific_enthalpy
    @specific_enthalpy.setter
    def specific_enthalpy(self):
        print(single_parameter_set_message)
    @specific_enthalpy.deleter
    def specific_enthalpy(self):
        raise AttributeError("Can't delete attribute") 
    # ---- Specific volume ----------------------------------------------------
    @property
    def specific_volume(self):
        """
        Specific volume of moist air, in m^3/kg.
        """
        return self._specific_enthalpy
    @specific_volume.setter
    def specific_volume(self):
        print(single_parameter_set_message)
    @specific_volume.deleter
    def specific_volume(self):
        raise AttributeError("Can't delete attribute")
    # ==== Couples of parameters to set the moist air physical properties =====
    # Dry temperature and relative humidity as entered parameters -------------
    @property
    def TR(self):
        """
        Dry temperature and relative humidity as entered properties of the moist
        air.
        """
        return self.temperature, self.relative_humidity
    @TR.setter
    def TR(self, values):
        # Set of the temperature value
        self._temperature = values[0]
        # Check if the entered value of relative humidity is between 0 and 1
        if (values[1] < 0) or (values[1] > 1):
            raise ValueError("Relative humidity must be between 0% and 100%!")
        else:
            self._relative_humidity = values[1]
        # Calculation of the corresponding specific humidity value
        self._specific_humidity = HAPropsSI('W', 'P', self.pressure*1e5,
                                            'T', self.temperature+273.15,
                                            'R', self.relative_humidity)
        # Wet bulb temperature
        self._wet_bulb_temperature = HAPropsSI('B', 'P', self.pressure*1e5,
                                               'T', self.temperature+273.15,
                                               'R', self.relative_humidity)\
                -273.15
        # Dew point temperature
        self._dew_point_temperature = HAPropsSI('D', 'P', self.pressure*1e5,
                                                'T', self.temperature+273.15,
                                                'R', self.relative_humidity)\
                -273.15
        # Specific enthalpy, in kJ/kg
        self._specific_enthalpy = HAPropsSI('H', 'P', self.pressure*1e5,
                                            'T', self.temperature+273.15,
                                            'R', self.relative_humidity)*1e-3
        # Specific volume, in m^3/kg
        self._specific_volume = HAPropsSI('V', 'P', self.pressure*1e5,
                                          'T', self.temperature+273.15,
                                          'R', self.relative_humidity)
    @property
    def RT(self):
        """
        Relative humidity and dry temperature as entered properties of the moist
        air.
        """
        return self.relative_humidity, self.temperature
    @RT.setter
    def RT(self, values):
        self.TR = values[::-1]
    # Dry and wet bulb temperature as entered parameters ----------------------
    @property
    def TB(self):
        """
        Dry and wet bulb temperatures as entered properties of the moist air.
        """
        return self.temperature, self.wet_bulb_temperature
    @TB.setter
    def TB(self, values):
        # Set of the temperature value
        self._temperature = values[0]
        # Check if the entered wet bulb temperature value is lower than the dry
        # one
        if (values[1] > values[0]):
            raise ValueError("Wet bulb temperature is lower than the dry one!")
        else:
            self._wet_bulb_temperature = values[1] 
        # The use of dry and wet bulb temperatures being usually quite slow with
        # CoolProp, the option fast_computation can be used here when a large
        # amount of temperature values is involved
        if not self.fast_computation:
            # Calculation of the corresponding relative humidity, using the
            # CoolProp package
            self._relative_humidity = HAPropsSI('R', 
                                                'P', self.pressure*1e5,
                                                'T', self.temperature+273.15,
                                                'B', self.wet_bulb_temperature\
                                                +273.15)
            # Calculation of the corresponding specific humidity value
            self._specific_humidity = HAPropsSI('W', 'P', self.pressure*1e5,
                                                'T', self.temperature+273.15,
                                                'B', self.wet_bulb_temperature\
                                                +273.15)
            # Dew point temperature
            self._dew_point_temperature = HAPropsSI('D', 'P', self.pressure*1e5,
                                                    'T',\
                                                    self.temperature+273.15,
                                                    'B',\
                                                    self.wet_bulb_temperature\
                                                    +273.15)-273.15
            # Specific enthalpy, in kJ/kg
            self._specific_enthalpy = HAPropsSI('H', 'P', self.pressure*1e5,
                                                'T', self.temperature+273.15,
                                                'B', self.wet_bulb_temperature\
                                                +273.15)*1e-3 
            # Specific volume, in m^3/kg
            self._specific_volume = HAPropsSI('V', 'P', self.pressure*1e5,
                                              'T', self.temperature+273.15,
                                              'B', self.wet_bulb_temperature\
                                              +273.15)
        else:
            # When fast computation is called, relative humidity is computed
            # in a simplified way, starting with the relative humidity 
            psatw = self.__equilibrium_vapor_pressure(self.wet_bulb_temperature)
            psatd = self.__equilibrium_vapor_pressure(self.temperature)
            beta = DRY_AIR_CP/(ALPHAW*WATER_LW)
            pvap = psatw-beta*self.pressure*1e5*\
                    (self.temperature-self.wet_bulb_temperature)
            self._relative_humidity = pvap/psatd
            # Specific humidity, see [1, page 6.10]
            self._specific_humidity = ALPHAW*pvap/(self.pressure*1e5-pvap)
            # Dew point temperature, see [1, page 6.10]
            lnpvap = np.log(pvap*1e-3)
            if self.temperature >= 0:
                self._dew_point_temperature = 6.54+14.526*lnpvap\
                        +0.7389*pow(lnpvap,2)+0.09486*pow(lnpvap,3)\
                        +0.4569*pow(pvap*1e-3,0.1984)
            else:
                self._dew_point_temperature = 6.09+12.608*lnpvap\
                        +0.4959*pow(lnpvap,2)
            # Specific enthalpy, in kJ/kg, see [1, page 6.9]
            self._specific_enthalpy = (DRY_AIR_CP*self.temperature\
                    +self.specific_humidity*(WATER_VAPOR_CP*self.temperature\
                                            +WATER_LW))*1e-3
            # Specific volume, in m^3/kg, using the ideal gas law
            self._specific_volume = WATER_VAPOR_R\
                    *(ALPHAW+self.specific_humidity)*(273.15+self.temperature)\
                    /((1+self.specific_humidity)*self.pressure*1e5)
    @property
    def BT(self):
        """
        Wet and dry bulb temperatures as entered properties of the moist air.
        """
        return self.wet_bulb_temperature, self.temperature
    @BT.setter
    def BT(self, values):
        self.TB = values[::-1]
    # Dry and dew point temperatures as entered properties --------------------
    @property
    def TD(self):
        """
        Dry and dew point temperatures as entered properties of the moist air.
        """
        return self.temperature, self.dew_point_temperature
    @TD.setter
    def TD(self, values):
        # Set of the temperature value
        self._temperature = values[0]
        # Check if the entered dew point temperature value is lower than the dry
        # one
        if (values[1] > values[0]):
            raise ValueError("Dew point temperature is lower than the dry one!")
        else:
            self._dew_point_temperature = values[1] 
        # Calculation of the corresponding relative humidity
        self._relative_humidity = HAPropsSI('R', 'P', self.pressure*1e5,
                                            'T', self.temperature+273.15,
                                            'D', self.dew_point_temperature\
                                            +273.15)
        # Calculation of the corresponding specific humidity value
        self._specific_humidity = HAPropsSI('W', 'P', self.pressure*1e5,
                                            'T', self.temperature+273.15,
                                            'D', self.dew_point_temperature\
                                            +273.15)
        # Dew point temperature
        self._wet_bulb_temperature = HAPropsSI('B', 'P', self.pressure*1e5,
                                               'T', self.temperature+273.15,
                                               'D', self.dew_point_temperature\
                                               +273.15)-273.15
        # Specific enthalpy, in kJ/kg
        self._specific_enthalpy = HAPropsSI('H', 'P', self.pressure*1e5,
                                            'T', self.temperature+273.15,
                                            'D', self.dew_point_temperature\
                                            +273.15)*1e-3
        # Specific volume, in m^3/kg
        self._specific_volume = HAPropsSI('V', 'P', self.pressure*1e5,
                                          'T', self.temperature+273.15,
                                          'D', self.dew_point_temperature\
                                          +273.15)
    @property
    def DT(self):
        """
        Dew point and dry bulb temperatures as entered properties of the moist
        air.
        """
        return self.dew_point_temperature, self.temperature
    @DT.setter
    def DT(self, values):
        self.TD = values[::-1]
    # Relative humidity and wet bulb temperature as entered properties --------
    @property
    def RB(self):
        """
        Relative humidity and wet bulb temperature as entered properties of the
        moist air.
        """
        return self.relative_humidity, self.wet_bulb_temperature
    @RB.setter
    def RB(self, values):
        # Set of the wet bulb temperature value
        self._wet_bulb_temperature = values[1]
        # Check if the entered value of relative humidity is between 0 and 1
        if (values[0] < 0) or (values[0] > 1):
            raise ValueError("Relative humidity must be between 0% and 100%!")
        else:
            self._relative_humidity = values[0]
        # Calculation of the corresponding dry bulb temperature
        self._temperature = HAPropsSI('T', 'P', self.pressure*1e5,
                                      'R', self.relative_humidity,
                                      'B', self.wet_bulb_temperature\
                                      +273.15)-273.15
        # Calculation of the corresponding specific humidity value
        self._specific_humidity = HAPropsSI('W', 'P', self.pressure*1e5,
                                            'R', self.relative_humidity,
                                            'B', self.wet_bulb_temperature\
                                            +273.15)
        # Dew point temperature
        self._dew_point_temperature = HAPropsSI('D', 'P', self.pressure*1e5,
                                                'R', self.relative_humidity,
                                                'B', self.wet_bulb_temperature\
                                                +273.15)-273.15
        # Specific enthalpy, in kJ/kg
        self._specific_enthalpy = HAPropsSI('H', 'P', self.pressure*1e5,
                                            'R', self.relative_humidity,
                                            'B', self.wet_bulb_temperature\
                                            +273.15)*1e-3 
        # Specific volume, in m^3/kg
        self._specific_volume = HAPropsSI('V', 'P', self.pressure*1e5,
                                          'R', self.relative_humidity,
                                          'B', self.wet_bulb_temperature\
                                          +273.15)
    @property
    def BR(self):
        """
        Wet bulb temperature and relative humidity as entered properties of the
        moist air.
        """
        return self.wet_bulb_temperature, self.relative_humidity
    @BR.setter
    def BR(self, values):
        self.RB = values[::-1]
    # Relative humidity and dew point temperature as entered properties --------
    @property
    def RD(self):
        """
        Relative humidity and dew point temperature as entered properties of the
        moist air.
        """
        return self.relative_humidity, self.dew_point_temperature
    @RD.setter
    def RD(self, values):
        # Set of the dew point temperature value
        self._dew_point_temperature = values[1]
        # Check if the entered value of relative humidity is between 0 and 1
        if (values[0] < 0) or (values[0] > 1):
            raise ValueError("Relative humidity must be between 0% and 100%!")
        else:
            self._relative_humidity = values[0]
        # Calculation of the corresponding dry bulb temperature
        self._temperature = HAPropsSI('T', 'P', self.pressure*1e5,
                                      'R', self.relative_humidity,
                                      'D', self.dew_point_temperature\
                                      +273.15)-273.15
        # Calculation of the corresponding specific humidity value
        self._specific_humidity = HAPropsSI('W', 'P', self.pressure*1e5,
                                            'R', self.relative_humidity,
                                            'D', self.dew_point_temperature\
                                            +273.15)
        # Wet bulb temperature
        self._wet_bulb_temperature = HAPropsSI('B', 'P', self.pressure*1e5,
                                               'R', self.relative_humidity,
                                               'D', self.dew_point_temperature\
                                               +273.15)-273.15
        # Specific enthalpy, in kJ/kg
        self._specific_enthalpy = HAPropsSI('H', 'P', self.pressure*1e5,
                                            'R', self.relative_humidity,
                                            'D', self.dew_point_temperature\
                                            +273.15)*1e-3
        # Specific volume, in m^3/kg
        self._specific_volume = HAPropsSI('V', 'P', self.pressure*1e5,
                                          'R', self.relative_humidity,
                                          'D', self.dew_point_temperature\
                                          +273.15)
    @property
    def DR(self):
        """
        Dew point temperature and relative humidity as entered properties of the
        moist air.
        """
        return self.dew_point_temperature, self.relative_humidity
    @DR.setter
    def DR(self, values):
        self.RD = values[::-1]
    # Wet bulb and dew point temperatures as entered properties ---------------
    @property
    def BD(self):
        """
        Wet bulb and dew point temperatures as entered properties of the moist
        air.
        """
        return self.wet_bulb_temperature, self.dew_point_temperature
    @BD.setter
    def BD(self, values):
        # Set of the wet bulb temperature value
        self._wet_bulb_temperature = values[0]
        # Set of the dew point temperature value
        self._dew_point_temperature = values[1]
        # Calculation of the corresponding dry bulb temperature
        self._temperature = HAPropsSI('T', 'P', self.pressure*1e5,
                                      'B', self.wet_bulb_temperature+273.15,
                                      'D', self.dew_point_temperature+273.15)\
                -273.15
        # Relative humidity     
        self._relative_humidity = HAPropsSI('R', 'P', self.pressure*1e5,
                                            'B', self.wet_bulb_temperature\
                                            +273.15,
                                            'D', self.dew_point_temperature\
                                            +273.15)
        # Calculation of the corresponding specific humidity value
        self._specific_humidity = HAPropsSI('W', 'P', self.pressure*1e5,
                                            'B', self.wet_bulb_temperature\
                                            +273.15,
                                            'D', self.dew_point_temperature\
                                            +273.15)
        # Specific enthalpy, in kJ/kg
        self._specific_enthalpy = HAPropsSI('H', 'P', self.pressure*1e5,
                                            'B', self.wet_bulb_temperature\
                                            +273.15,
                                            'D', self.dew_point_temperature\
                                            +273.15)*1e-3
        # Specific volume, in m^3/kg
        self._specific_volume = HAPropsSI('V', 'P', self.pressure*1e5,
                                          'B', self.wet_bulb_temperature+273.15,
                                          'D', self.dew_point_temperature\
                                          +273.15)
    @property
    def DB(self):
        """
        Dew point and wet bulb temperatures humidity as entered properties of
        the moist air.
        """
        return self.dew_point_temperature, self.wet_bulb_temperature
    @DB.setter
    def DB(self, values):
        self.BD = values[::-1]
    # Dry temperature and specific humidity as entered properties -------------
    @property
    def TW(self):
        """
        Dry temperature and specific humidity as entered properties of the moist
        air.
        """
        return self.temperature, self.specific_humidity
    @TW.setter
    def TW(self, values):
        # Set of the dry temperature value
        self._temperature = values[0] 
        # Check if the entered value of specific humidity makes sense, so if it
        # is between 0 and its maximum value when saturation is reached at the
        # same dry temperature
        wmax = HAPropsSI('W', 'P', self.pressure*1e5, 'R', 1.0,\
                         'T', values[0]+273.15)
        if (values[1] < 0):
            raise ValueError("Specific humidity is positive")
        elif (values[1] > wmax):
            raise ValueError("Specific humidity is higher than at saturation")
        else:
            # Set of the specific humidity value
            self._specific_humidity = values[1]
        # Relative humidity     
        self._relative_humidity = HAPropsSI('R', 'P', self.pressure*1e5,
                                            'T', self.temperature+273.15,
                                            'W', self.specific_humidity)
        # Calculation of the corresponding wet bulb temperature 
        self._wet_bulb_temperature = HAPropsSI('B', 'P', self.pressure*1e5,
                                               'T', self.temperature+273.15,
                                               'W',
                                               self.specific_humidity)-273.15
        # Calculation of the corresponding dew point temperature 
        self._dew_point_temperature = HAPropsSI('D', 'P', self.pressure*1e5,
                                                'T', self.temperature+273.15,
                                                'W', self.specific_humidity)\
                -273.15
        # Specific enthalpy, in kJ/kg
        self._specific_enthalpy = HAPropsSI('H', 'P', self.pressure*1e5,
                                            'T', self.temperature+273.15,
                                            'W', self.specific_humidity)*1e-3
        # Specific volume, in m^3/kg
        self._specific_volume = HAPropsSI('V', 'P', self.pressure*1e5,
                                          'T', self.temperature+273.15,
                                          'W', self.specific_humidity)
    @property
    def WT(self):
        """
        Specific humidity and dry temperature as entered properties of the moist
        air.
        """
        return self.specific_humidity, self.temperature
    @WT.setter
    def WT(self, values):
        self.TW = values[::-1]
    # Relative and specific humidities as entered parameters ------------------
    @property
    def RW(self):
        """
        Relative and specific humidities as entered properties of the moist
        air.
        """
        return self.relative_humidity, self.specific_humidity
    @RW.setter
    def RW(self, values):
        # Check if the entered value of relative humidity is between 0 and 1
        if (values[0] < 0) or (values[0] > 1):
            raise ValueError("Relative humidity must be between 0% and 100%!")
        else:
            self._relative_humidity = values[0]
        # Check if the entered value of specific humidity makes sense
        if (values[1] < 0):
            raise ValueError("Specific humidity is positive")
        else:
            # Set of the specific humidity value
            self._specific_humidity = values[1]
        # Dry temperature
        self._temperature = HAPropsSI('T', 'P', self.pressure*1e5,
                                      'R', self.relative_humidity,
                                      'W', self.specific_humidity)-273.15
        # Calculation of the corresponding wet bulb temperature 
        self._wet_bulb_temperature = HAPropsSI('B', 'P', self.pressure*1e5,
                                               'R', self.relative_humidity,
                                               'W', self.specific_humidity)\
                -273.15
        # Calculation of the corresponding dew point temperature 
        self._dew_point_temperature = HAPropsSI('D', 'P', self.pressure*1e5,
                                                'R', self.relative_humidity,
                                                'W', self.specific_humidity)\
                -273.15
        # Specific enthalpy, in kJ/kg
        self._specific_enthalpy = HAPropsSI('H', 'P', self.pressure*1e5,
                                            'R', self.relative_humidity,
                                            'W', self.specific_humidity)*1e-3
        # Specific volume, in m^3/kg
        self._specific_volume = HAPropsSI('V', 'P', self.pressure*1e5,
                                          'R', self.relative_humidity,
                                          'W', self.specific_humidity)
    @property
    def WR(self):
        """
        Specific and relative humidities as entered properties of the moist
        air.
        """
        return self.specific_humidity, self.temperature
    @WR.setter
    def WR(self, values):
        self.RW = values[::-1] 
    # Wet bulb temperature and specific humidity as entered properties -------------
    @property
    def BW(self):
        """
        Wet bulb temperature and specific humidity as entered properties of the
        moist air.
        """
        return self.wet_bulb_temperature, self.specific_humidity
    @BW.setter
    def BW(self, values):
        # Set of the wet bulb temperature value
        self._wet_bulb_temperature = values[0] 
        # Check if the entered value of specific humidity makes sense, so if it
        # is between 0 and its maximum value when saturation is reached at the
        # same dry temperature
        wmax = HAPropsSI('W', 'P', self.pressure*1e5, 'R', 1.0,\
                         'B', values[0]+273.15)
        if (values[1] < 0):
            raise ValueError("Specific humidity is positive")
        elif (values[1] > wmax):
            raise ValueError("Specific humidity is higher than at saturation")
        else:
            # Set of the specific humidity value
            self._specific_humidity = values[1]
        # Dry temperature 
        self._temperature = HAPropsSI('T', 'P', self.pressure*1e5,
                                      'B', self.wet_bulb_temperature+273.15,
                                      'W', self.specific_humidity)-273.15
        # Relative humidity     
        self._relative_humidity = HAPropsSI('R', 'P', self.pressure*1e5,
                                            'B', self.wet_bulb_temperature\
                                            +273.15,
                                            'W', self.specific_humidity)
        # Calculation of the corresponding dew point temperature 
        self._dew_point_temperature = HAPropsSI('D', 'P', self.pressure*1e5,
                                                'B', self.wet_bulb_temperature\
                                                +273.15,
                                                'W', self.specific_humidity)\
                -273.15
        # Specific enthalpy, in kJ/kg
        self._specific_enthalpy = HAPropsSI('H', 'P', self.pressure*1e5,
                                            'B', self.wet_bulb_temperature\
                                            +273.15,
                                            'W', self.specific_humidity)*1e-3
        # Specific volume, in m^3/kg
        self._specific_volume = HAPropsSI('V', 'P', self.pressure*1e5,
                                          'B', self.wet_bulb_temperature+273.15,
                                          'W', self.specific_humidity)
    @property
    def WB(self):
        """
        Specific humidity and wet bulb temperature as entered properties of the
        moist air.
        """
        return self.specific_humidity, self.wet_bulb_temperature
    @WB.setter
    def WB(self, values):
        self.BW = values[::-1]
    # Dew point temperature and specific humidity as entered properties -------
    @property
    def DW(self):
        """
        Dew point temperature and specific humidity as entered properties(self.specific_humidity*WATER_VAPOR_R\
                *np.log(self.specific_humidity/ref_state.specific_humidity)\
                -(DRY_AIR_R+self.specific_humidity of the
        moist air.
        """
        return self.dew_point_temperature, self.specific_humidity
    @DW.setter
    def DW(self, values):
        # Set of the dew point temperature value
        self._dew_point_temperature = values[0] 
        # Check if the entered value of specific humidity makes sense, so if it
        # is between 0 and its maximum value when saturation is reached at the
        # same dry temperature
        wmax = HAPropsSI('W', 'P', self.pressure*1e5, 'R', 1.0,\
                         'D', values[0]+273.15)
        if (values[1] < 0):
            raise ValueError("Specific humidity is positive")
        elif (values[1] > wmax):
            raise ValueError("Specific humidity is higher than at saturation")
        else:
            # Set of the specific humidity value
            self._specific_humidity = values[1]
        # Dry temperature 
        self._temperature = HAPropsSI('T', 'P', self.pressure*1e5,
                                      'D', self.dew_point_temperature+273.15,
                                      'W', self.specific_humidity)-273.15
        # Relative humidity     
        self._relative_humidity = HAPropsSI('R', 'P', self.pressure*1e5,
                                            'D', self.dew_point_temperature\
                                            +273.15,
                                            'W', self.specific_humidity)
        # Calculation of the corresponding wet bulb temperature 
        self._wet_bulb_temperature = HAPropsSI('B', 'P', self.pressure*1e5,
                                               'D', self.dew_point_temperature\
                                               +273.15,
                                               'W', self.specific_humidity)\
                -273.15
        # Specific enthalpy, in kJ/kg
        self._specific_enthalpy = HAPropsSI('H', 'P', self.pressure*1e5,
                                            'D', self.dew_point_temperature\
                                            +273.15,
                                            'W', self.specific_humidity)*1e-3
        # Specific volume, in m^3/kg
        self._specific_volume = HAPropsSI('V', 'P', self.pressure*1e5,
                                          'D', self.dew_point_temperature\
                                          +273.15,
                                          'W', self.specific_humidity)
    @property
    def WD(self):
        """
        Specific humidity and dew point temperature as entered properties of the
        moist air.
        """
        return self.specific_humidity, self.dew_point_temperature
    @WD.setter
    def WD(self, values):
        self.DW = values[::-1]
    # ==== Usual class methods ================================================
    def specific_thermal_exergy(self, ref_state):
        """
        Specific thermal exergy of moist air, in kJ/kg, regarding to a reference
        state represented by its corresponding MoistAir object.
        """
        # Check if the entered reference state is a MoistAir object
        if type(ref_state) != MoistAir:
            raise TypeError("Reference state must be a MoisAir object!")
        # Calculation of the Carnot factor
        Tml0 = self.__logmeantemperature(self.temperature,\
                                         ref_state.temperature)
        CarnotFactor = 1-(ref_state.temperature+273.15)/Tml0
        # And specific thermal exergy
        return 1e-3*(DRY_AIR_CP+self.specific_humidity*WATER_VAPOR_CP)\
                *(self.temperature-ref_state.temperature)*CarnotFactor
    def specific_mechanical_exergy(self, ref_state):
        """
        Specific mechanical exergy of moist air, in kJ/kg, regarding to a
        reference state represented by its corresponding MoistAir object.
        """
        # Check if the entered reference state is a MoistAir object
        if type(ref_state) != MoistAir:
            raise TypeError("Reference state must be a MoisAir object!")
        return (DRY_AIR_R+self.specific_humidity*WATER_VAPOR_R)\
                *(ref_state.temperature+273.15)\
                *np.log(self.pressure/ref_state.pressure)*1e-3
    def specific_chemical_exergy(self, ref_state):
        """
        Specific chemical exergy of moist air, in kJ/kg, regarding to a
        reference state represented by its corresponding MoistAir object.
        """
        # Check if the entered reference state is a MoistAir object
        if type(ref_state) != MoistAir:
            raise TypeError("Reference state must be a MoisAir object!")
        # Specific humidities of each moist air
        w , w0 = self.specific_humidity, ref_state.specific_humidity
        T0 = ref_state.temperature+273.15
        return T0*(w*WATER_VAPOR_R*np.log(w/w0)-(DRY_AIR_R+w*WATER_VAPOR_R)\
                   *np.log((w+ALPHAW)/(w0+ALPHAW)))*1e-3
    def specific_exergy(self, ref_state):
        """
        Specific exergy of moist air, in kJ/kg, regarding to a reference state
        represented by its corresponding MoistAir object.
        """
        # Check if the entered reference state is a MoistAir object
        if type(ref_state) != MoistAir:
            raise TypeError("Reference state must be a MoisAir object!")
        return self.specific_thermal_exergy(ref_state)\
                +self.specific_mechanical_exergy(ref_state)\
                +self.specific_chemical_exergy(ref_state)


if __name__ == '__main__':
    inside = MoistAir()
    inside.TW = 20.0, 7e-3
    outside = MoistAir()
    outside.TW = 0.0, 3e-3

# ==== References =============================================================
#
# [1] “ASHRAE Handbook : Fundamentals” (2001). In : American Society of Heating,
# Refrigerating and Air-Conditioning Engineers. Chap. 6 : Psychrometrics.
