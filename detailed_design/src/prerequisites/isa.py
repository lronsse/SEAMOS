"""
Simple ISA-Calculator (written by Lilly Zwikker)
"""

#  --- Constants ---
gravity_sea_level = 9.80665  # m/s2
gas_constant = 287.0  # J/kgK

#  --- Sea level values ---

temperature_sea_level = 288.15  # K
pressure_sea_level = 101325  # Pa
density_sea_level = 1.225  # rho

# Layer LUT
layer_boundaries = [0, 11000, 20000, 32000, 47000, 51000, 71000, 86000]
layer_gradients = [-0.0065, 0, 0.0010, 0.0028, 0, -0.0028, -0.0020]


class HeightError(Exception):
    """
    Error for invalid heights
    """
    pass


class MenuError(Exception):
    """
    Error for invalid menu entries
    """
    pass


def temperature_kelvin_to_celsius(temp):
    """
    :param temp: Temperature in Kelvin
    :return: Temperature in Celsius
    """
    return temp - 273.15


def temperature_celsius_to_kelvin(temp):
    """
    :param temp: Temperature in Celsius
    :return: Temperature in Kelvin
    """
    return temp + 273.15


def ceiling(value: float | int, ceiling_value: float | int) -> float | int:
    """
    :param value: Value that cannot exceed ceiling
    :param ceiling_value: Ceiling value that can not be exceeded
    :return: limited value
    """
    if value <= ceiling_value:
        return value
    if value > ceiling_value:
        return ceiling_value


def calculate_temperature(temperature_start: float, gradient: float, height_start: float, height_end: float) -> float:
    """
    calculate new temperature based on constant temperature gradient

    :param temperature_start: temperature in Celsius
    :param gradient: constant temperature gradient in Celsius/meter
    :param height_start: start height in meter
    :param height_end: target height in meter
    :return: temperature in Celsius
    """
    return temperature_start + gradient * (height_end - height_start)


def calculate_new_pressure(pressure_start: float, temperature_start: float, temperature_end: float, gradient: float,
                           height_start: float, height_end: float) -> float:
    """

    :param pressure_start: pressure in Pa
    :param temperature_start: temperature in Celsius
    :param temperature_end: temperature in Celsius
    :param gradient: constant temperature gradient in Celsius/meter
    :param height_start: height in meter
    :param height_end: height in meter
    :return: pressure in Pa
    """
    if gradient != 0:
        return pressure_start * (temperature_end / temperature_start) ** (-gravity_sea_level / (gradient * gas_constant))
    if gradient == 0:
        return pressure_start * 2.71828 ** (-(gravity_sea_level / (gas_constant * temperature_start))
                                            * (height_end - height_start))


def calculate_density(pressure_local: float, temperature_local: float) -> float:
    """
    calculate the density based on the gas constant, temperature and pressure

    :param pressure_local: pressure in Pa
    :param temperature_local: temperature in Celsius
    :return: density in kg/m3
    """
    return pressure_local / (gas_constant * temperature_local)


def correct_units(input_value: float, unit: str) -> float:
    """
    convert from meter, feet and FL to meter

    :param input_value: the value to convert
    :param unit: 'meter', 'feet', 'FL'
    :return: the input value converted to the right units
    """
    input_value = int(input_value)
    if unit == 'meter':
        return input_value
    if unit == 'feet':
        return input_value / 3.281
    if unit == 'FL':
        return (input_value * 100) / 3.281
    else:
        raise MenuError


def calculate_isa(height_input: float) -> (float, float, float):
    """
    calculate the ISA values at specific height by iterating over layers

    :param height_input: height in meters
    :return: pressure in Pa, temperature in Celsius, density in kg/m3
    """
    current_temperature = temperature_sea_level
    current_pressure = pressure_sea_level
    if height_input <= 0 or height_input > 86000:
        raise HeightError
    for layer_i in range(1, len(layer_boundaries)):
        if height_input > layer_boundaries[layer_i - 1]:
            local_height = ceiling(height_input, layer_boundaries[layer_i])
            old_temperature = current_temperature
            current_temperature = calculate_temperature(old_temperature, layer_gradients[layer_i - 1],
                                                        layer_boundaries[layer_i - 1], local_height)
            current_pressure = calculate_new_pressure(current_pressure, old_temperature, current_temperature,
                                                      layer_gradients[layer_i - 1], layer_boundaries[layer_i - 1],
                                                      local_height)

    current_density = calculate_density(current_pressure, current_temperature)

    rounded_pressure = round(current_pressure, 3)
    rounded_temperature = round(current_temperature, 3)
    rounded_density = round(current_density, 5)

    return rounded_pressure, rounded_temperature, rounded_density


def main() -> None:
    """
    UI wrapper for ISA Calculator
    """
    try:
        print("*** ISA Calculator ***\n")
        print("1. Calculate ISA for altitude in meters\n2. Calculate ISA for altitude in feet\n3. Calculate ISA for "
              "altitude in FL\n")
        unit_choice = input("Enter your choice: ")
        if unit_choice not in ['meter', 'feet', 'FL']:
            raise MenuError("Unit is entered incorrectly or not supported")
        height = int(input("Height: "))
        height = correct_units(height, unit_choice)
        pressure, temperature, density = calculate_isa(height)
        print(" Pressure:", pressure, "\n Temperature: ", temperature, "\n Density: ", density)
    except HeightError:
        print("ERROR: Height is too high or too low")
    except ValueError:
        print("ERROR: You entered a wrong type of value")
    except MenuError:
        print("ERROR: That is not a valid option")


if __name__ == '__main__':
    main()

"""
Tests for ISA Calculator
"""
import pytest


def test_unit_converters():
    assert temperature_kelvin_to_celsius(0) == -273.15
    assert temperature_celsius_to_kelvin(0) == 273.15
    assert correct_units(10, 'meter') == 10
    assert correct_units(10, 'feet') == 10 / 3.281
    assert correct_units(10, 'feet') != 10
    assert correct_units(10, 'FL') != 10
    assert correct_units(10, 'FL') == 1000 / 3.281


@pytest.mark.parametrize("test_value,test_ceiling,expected_result", [(10, 20, 10), (30, 20, 20), (20, 20, 20), (-5, 20, -5),
                                                                     (-10, -20, -20), (20, 10, 10)])
def test_ceiling_function(test_value, test_ceiling, expected_result):
    assert ceiling(test_value, test_ceiling) == expected_result


@pytest.mark.parametrize("test_height,expected_result", [(1, (101312.985, 288.143, 1.22511)),
                                                         (85999, (0.302, 184.652, 1e-05))])
def test_calculate_isa(test_height, expected_result):
    assert calculate_isa(test_height) == expected_result


@pytest.mark.parametrize("test_height", [-1, 0, 86001])
def test_calculate_isa_height(test_height):
    with pytest.raises(HeightError):
        calculate_isa(test_height)
