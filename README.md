# PsychroPy: Psychrometrics with Python

## Presentation

**PsychroPy** is a python package dedicated to
[psychrometrics](https://en.wikipedia.org/wiki/Psychrometrics) calculations, so
dealing with _humid_ or _moist air_. This package is (so far) mostly based on
the [CoolProp](http://www.coolprop.org/) one, and more specifically on its
[humid air](http://www.coolprop.org/fluid_properties/HumidAir.html) function
called `HAPropsSI`. The main brick of this package is for the moment a class
called `MoistAir`, which contain the complete physical state of a given humid
air, so its **temperature** and its
[**humidity**](https://en.wikipedia.org/wiki/Humidity), considered for example
at a specific location within an [HVAC](https://en.wikipedia.org/wiki/HVAC)
system. The humidity of air can be entered using either:
- its _specific humidity_ (aka _moisture content_ or _mixing ratio_), noted `W`
  and defined as the ratio of the mass of _water vapour_ on the one of **dry
air** it is mixed with, so mathematically: $$W = \frac{m_w}{m_{da}}$$
- its [_relative humidity_](https://en.wikipedia.org/wiki/Relative_humidity), noted `R`.
- its [_wet bulb temperature_](https://en.wikipedia.org/wiki/Wet-bulb_temperature), noted `B`.
- its [_dew point temperature_](https://en.wikipedia.org/wiki/Dew_point), noted `D`.
From such parameters, all the other ones, as for instance the air
_specific enthalpy_ or its _specific volume_, can be directly obtained. Two
parameters, as e.g. _temperature_ and _relative humidity_ (`TR`) or _dry_ and
_wet bulb temperature_ (`TB`), are **required** to unambiguously define the
state of any humid air.

## Example

Consider for example an hypothetical HVAC system that extracts air from the
outside and, after treatment, inject it into a given building in order to reach
a humid air state corresponding to a temperature of 20°C and a relative
humidity of 50%. To manipulate such an air, we can instance an object from the
class `MoistAir`, called for example `inside`, which will contain what we know
about this air: 
```python
from psychrometry import psychropy as psp
inside =  psp.MoistAir()
# Known parameters of outside air
inside.TR = 20.0 + 273.15, 0.5
# Requested specific humidity and enthalpy of the inside
Wi, Hi = inside.specific_humidity, inside.specific_enthalpy
```
The results are a specific humidity `Wi = 7.29e-3` and a specific enthalpy `Hi = 38.62` (in kJ/kg). Let us emphasize that:
- Temperature is always expressed in [Kelvin](https://en.wikipedia.org/wiki/Kelvin).
- Parameters of the air are entered using a notation close the one used by [CANTERA](https://cantera.org/), and mentioning explicitely the two concerned parameters, here `inside.TR` for temperature and relative humidity. The same parameters may be entered in the reverse order using `inside.RT`, or using different ways to quantity humidity: `inside.TB`, `inside.DT`, `inside.TW`, and so on.
- Specific humidity is expressed in kg of water vapour per kilogram of dry air, so dimensionless. 
- Specific enthalpy is expressed in kJ/kg or dry air.
If the HVAC system is equipped with two temperature sensors, a usual one (that gives `T`) and a [wet bulb](https://en.wikipedia.org/wiki/Wet-bulb_temperature) one (that gives `B`), both measuring the state of the outside air. We can have for example (in winter):
```python
outside = psp.MoistAir()
# Outside pressure of 1 bar
outside.pressure = 1.0 * 1e5
# Measured dry temperature of 5°C and Wet-Bulb one of 3°C
outside.TB = 5.0 + 273.15, 3.0 + 273.15
# Outside specific humidity and enthalpy 
W0, H0 = outside.specific_humidity, outside.specific_enthalpy
```
If the dry air mass flow rate of ventilation is known, and noted for example `ma`, the humidifying water mass flow rate to inject into this air stream is given by:
```python
# Humidification mass flow rate, with the same unit of the dry air one 
humidification_flow_rate = ma*(Wi-W0)
```
In the very same way, the heat rate with have to supply to the same air stream is given by:
```python
# Heating heat rate, expressed in kW if the dry air mass flow rate is in kg/s
heat_rate = ma*(Hi-H0)
```

## Exergy calculations ##

[**Exergy**](https://en.wikipedia.org/wiki/Exergy) of humid air, regarding to a
_reference_ or _dead state_ can also be calculated by **PsychroPy**, using:
```python
ex = inside.specific_exergy(ref_state=inside)
```
The used dead state has to be represented by an `MoistAir` object.

## The _fast computation_ option

Although very accurate and efficient, the [CoolProp](http://www.coolprop.org/) package involves sometimes quite long calculations, for example when the air parameters are entered using the [wet bulb](https://en.wikipedia.org/wiki/Wet-bulb_temperature) or the [dew_point](https://en.wikipedia.org/wiki/Dew_point) temperature. When large amounts of temperature measurements are involved, it is possible to switch the `psychropy` package into a **fast computation** mode, using the `fast_computation` boolean attribute, as:
```python
outside.fast_computation = True
```
The [CoolProp](http://www.coolprop.org/) package is then not used anymore for calculations, which are led using the common explicit equations of psychrometry, as presented in the well known [ASHRAE](https://www.ashrae.org/) [handbook](https://www.ashrae.org/technical-resources/ashrae-handbook).

Any question, remark or advice of improvement is welcome.

