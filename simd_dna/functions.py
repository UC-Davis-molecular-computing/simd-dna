def convert_hex_to_rgb(hex_rgb: str) -> str:
    """Returns an RGB string representation of a hexadecimal color code

    :param hex_rgb: A 6-digit hex color code
    :return: An RGB string representation of the provided hex color code
    """
    if hex_rgb[0] != '#':
        hex_rgb = '#' + hex_rgb

    red = int(hex_rgb[1:3], 16)
    green = int(hex_rgb[3:5], 16)
    blue = int(hex_rgb[5:7], 16)
    return 'rgb(%d, %d, %d)' % (red, green, blue)


def convert_rgb_to_hex(r: float, g: float, b: float) -> str:
    """Returns a hexadecimal color code representation of a set of RGB values

    :param r: A float representing the red component
    :param g: A float representing the green component
    :param b: A float representing the blue component
    :return: A 6-digit hex color code representing the provided RGB values
    """
    red = int(r * 255)
    green = int(g * 255)
    blue = int(b * 255)
    return '#{:02x}{:02x}{:02x}'.format(red, green, blue)
