def convert_hex_to_rgb(hex_rgb):
    red = int(hex_rgb[1:3], 16)
    green = int(hex_rgb[3:5], 16)
    blue = int(hex_rgb[5:7], 16)
    return 'rgb(%d, %d, %d)' % (red, green, blue)


def convert_rgb_to_hex(r, g, b):
    red = int(r * 255)
    green = int(g * 255)
    blue = int(b * 255)
    return '#{:02x}{:02x}{:02x}'.format(red, green, blue)
