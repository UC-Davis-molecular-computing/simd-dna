import copy
import sys
import json
import svgwrite
import re
from json import JSONEncoder

program_loop = True
step_by_step_simulation = False
keep_results = False
show_unused_instruction_strands = False
cell_types = {}
strand_types = {}
registers = {}
instructions = []


class ObjectEncoder(JSONEncoder):
    def default(self, o):
        return o.__dict__


class Strand:

    def __init__(self, domains, is_complementary, color='#000000'):
        self.domains = domains
        self.is_complementary = is_complementary
        if re.match('^#?[A-Fa-f0-9]{6}$', color) is not None:
            self.color = color
            if not self.color.startswith('#'):
                self.color = '#' + self.color
        else:
            self.color = '#000000'

    @staticmethod
    def decode_json(domains, is_complementary, color='#000000'):
        self = Strand(domains, is_complementary, color)
        return self


class Cell:

    def __init__(self, domains):
        self.domains = domains

    @staticmethod
    def decode_json(domains):
        self = Cell(domains)
        return self


class Register:

    def __init__(self):
        self.cells = []
        self.coverings = []
        self._dwg = None
        self._svg_domain_length = 5     # millimeters
        self._svg_left_offset = 5
        self._svg_vertical_offset = 55
        self._svg_cell_height = 40
        self._total_domains = 0

    def add_cell(self, cell_name):
        self.cells.append(cell_name)
        self._total_domains += cell_types[cell_name]["domains"]

    def get_cell_at_domain_index(self, domain_index):
        total_domains = 0
        for cell_name in self.cells:
            cell = cell_types[cell_name]
            if total_domains <= domain_index < total_domains + len(cell.domains):
                return cell, domain_index - total_domains
            total_domains += len(cell.domains)

        return None, 0

    def get_coverings_at_domain_index(self, domain_index, include_orthogonal=False, strand_set=None):
        coverings = []
        cell, offset = self.get_cell_at_domain_index(domain_index)
        if cell is None:
            return coverings
        domain_label = cell.domains[offset]

        if include_orthogonal:
            orthogonal_coverings = []

        if strand_set is None:
            strand_set = self.coverings
        for covering in strand_set:
            start_index = covering['start_index']
            strand = strand_types[covering['strand_name']]
            if start_index <= domain_index < start_index + len(strand.domains):
                if domain_label == strand.domains[domain_index - start_index]:
                    coverings.append(covering)
                elif include_orthogonal:
                    orthogonal_coverings.append(covering)

        if include_orthogonal:
            return coverings, orthogonal_coverings
        else:
            return coverings

    def attempt_attachment(self, domain_index, strand_type, unattached_matches=None):
        if domain_index < 0:
            total_domains = 0
            for cell_name in self.cells:
                total_domains += len(cell_types[cell_name].domains)

            domain_index += total_domains
        strand = strand_types[strand_type]

        if strand.is_complementary:
            displaced_strands = []
            displacing_strands = []
            for covering in self.coverings:
                top_strand = strand_types[covering['strand_name']]
                is_match = True
                for i in range(len(top_strand.domains)):
                    if strand.domains[i] != top_strand.domains[i]:
                        is_match = False
                        break

                if is_match:
                    strand_start = covering['start_index']
                    strand_end = strand_start + len(top_strand.domains)
                    for i in range(strand_start, strand_end):
                        coverings_at_domain = self.get_coverings_at_domain_index(i)
                        # must have at least one insecure domain
                        if len(coverings_at_domain) > 1 or covering not in coverings_at_domain:
                            displaced_strands.append(covering)
                            displacing_strands.append({'start_index': strand_start, 'strand_name': strand_type})
                            break

            if len(displaced_strands) > 0:
                self.coverings = [x for x in self.coverings if x not in displaced_strands]
                return displacing_strands
            elif unattached_matches is not None:
                matchings = 0
                for i in range(len(strand.domains)):
                    cell, offset = self.get_cell_at_domain_index(domain_index + i)
                    if cell is not None and cell.domains[offset] == strand.domains[i]:
                        matchings += 1

                if matchings >= 1:
                    new_strand = {'start_index': domain_index, 'strand_name': strand_type}
                    if new_strand not in unattached_matches:
                        unattached_matches.append(new_strand)
                        unattached_matches.sort(key=lambda x: x['start_index'])
        else:
            has_open_toehold = False
            for i in range(len(strand.domains)):
                coverings = self.get_coverings_at_domain_index(domain_index + i)
                if len(coverings) == 0:
                    cell, offset = self.get_cell_at_domain_index(domain_index + i)
                    if cell is not None and cell.domains[offset] == strand.domains[i]:
                        has_open_toehold = True

            if has_open_toehold or unattached_matches is not None:
                # must have at least two matching domains to attach
                matchings = 0
                for i in range(len(strand.domains)):
                    cell, offset = self.get_cell_at_domain_index(domain_index + i)
                    if cell is not None and cell.domains[offset] == strand.domains[i]:
                        matchings += 1

                if matchings >= 2:
                    new_covering = {'start_index': domain_index, 'strand_name': strand_type}
                    if has_open_toehold:
                        self.coverings.append(new_covering)
                        self.coverings.sort(key=lambda x: x['start_index'])
                        return [new_covering]
                    else:
                        if new_covering not in unattached_matches:
                            unattached_matches.append(new_covering)
                            unattached_matches.sort(key=lambda x: x['start_index'])
                        return None

        return None

    def displace_strands(self, coverings=[]):
        displaced_strands = []
        coverings = [x for x in self.coverings if x not in coverings]

        for covering in coverings:
            strand_start = covering['start_index']
            strand = strand_types[covering['strand_name']]
            strand_end = strand_start + len(strand.domains)
            insecure_domains = 0
            for i in range(strand_start, strand_end):
                coverings_at_domain = self.get_coverings_at_domain_index(i)
                if len(coverings_at_domain) > 1 or covering not in coverings_at_domain:
                    insecure_domains += 1

            if insecure_domains >= strand_end - strand_start - 1:
                displaced_strands.append(covering)

        if len(displaced_strands) > 0:
            self.coverings = [x for x in self.coverings if x not in displaced_strands]

        return displaced_strands

    def print(self, new_strands=None, unused_strands=None):
        if unused_strands is not None and len(unused_strands) > 0:
            self._print_floating_strands(unused_strands)

        if new_strands is not None and len(new_strands) > 0:
            self._print_floating_strands(new_strands)
        elif unused_strands is not None and len(unused_strands) > 0:
            self._print_empty_layer('-')
        else:
            self._print_empty_layer()

        previous_domains = 0
        print('|', end='')
        for cell_name in self.cells:
            cell = cell_types[cell_name]
            for i in range(len(cell.domains)):
                domain_coverings, orthogonal_coverings = self.get_coverings_at_domain_index(previous_domains + i,
                                                                                            include_orthogonal=True)
                if len(orthogonal_coverings) >= 1:
                    point_right = True
                    for covering in orthogonal_coverings:
                        if covering['start_index'] == previous_domains + i:
                            point_right = False
                            break

                    if point_right:
                        print('/', end='')
                    else:
                        print('\\', end='')
                else:
                    print(' ', end='')

            print('|', end='')
            previous_domains += len(cell.domains)

        if len(self.coverings) >= 1:
            last_covering = self.coverings[-1]
            strand = strand_types[last_covering['strand_name']]
            for _ in range(previous_domains, last_covering['start_index'] + len(strand.domains)):
                print('/', end='')

        print()

        previous_domains = 0
        print('|', end='')
        for cell_name in self.cells:
            cell = cell_types[cell_name]
            for i in range(len(cell.domains)):
                coverings = self.get_coverings_at_domain_index(previous_domains + i)
                if len(coverings) == 0:
                    print('□', end='')
                elif len(coverings) == 1:
                    strand = strand_types[coverings[0]['strand_name']]
                    index = previous_domains + i - coverings[0]['start_index']
                    if index < len(strand.domains) - 1:
                        print('=', end='')
                    else:
                        print('>', end='')
                else:
                    print('x', end='')

            print('|', end='')
            previous_domains += len(cell.domains)

        print()

    def _print_floating_strands(self, strand_set):
        previous_domains = 0
        print('|', end='')
        for cell_name in self.cells:
            cell = cell_types[cell_name]
            for i in range(len(cell.domains)):
                domain_coverings, orthogonal_coverings = self.get_coverings_at_domain_index(previous_domains + i,
                                                                                            include_orthogonal=True,
                                                                                            strand_set
                                                                                            =strand_set)
                if len(orthogonal_coverings) >= 1:
                    point_right = True
                    for covering in orthogonal_coverings:
                        if covering['start_index'] == previous_domains + i:
                            point_right = False
                            break

                    if point_right:
                        print('/', end='')
                    else:
                        print('\\', end='')
                else:
                    print(' ', end='')

            print('|', end='')
            previous_domains += len(cell.domains)

        if len(strand_set) >= 1:
            last_covering = strand_set[-1]
            strand = strand_types[last_covering['strand_name']]
            for _ in range(previous_domains, last_covering['start_index'] + len(strand.domains)):
                print('/', end='')

        print()

        previous_domains = 0
        print('|', end='')
        for cell_name in self.cells:
            cell = cell_types[cell_name]
            for i in range(len(cell.domains)):
                coverings = self.get_coverings_at_domain_index(previous_domains + i, strand_set=strand_set)
                if len(coverings) == 0:
                    print(' ', end='')
                elif len(coverings) == 1:
                    strand = strand_types[coverings[0]['strand_name']]
                    index = previous_domains + i - coverings[0]['start_index']
                    if strand.is_complementary:
                        if index == 0:
                            print('<', end='')
                        else:
                            print('=', end='')
                    else:
                        if index < len(strand.domains) - 1:
                            print('=', end='')
                        else:
                            print('>', end='')
                else:
                    print('x', end='')

            print('|', end='')
            previous_domains += len(cell.domains)

        print()

    def _print_empty_layer(self, blank_char=' '):
        previous_domains = 0
        print('|', end='')
        for cell_name in self.cells:
            cell = cell_types[cell_name]
            for i in range(len(cell.domains)):
                print(blank_char, end='')

            print('|', end='')
            previous_domains += len(cell.domains)

        print()

    def sanitize_unattached_strands(self, unattached_strands, new_strands):
        sanitized_strands = []
        for strand in reversed(unattached_strands):
            add_strand = True
            for new_strand in new_strands:
                if self.strands_intersect(strand, new_strand):
                    add_strand = False
                    break

            if add_strand:
                for other_strand in sanitized_strands:
                    if self.strands_intersect(strand, other_strand):
                        add_strand = False
                        break

            if add_strand:
                sanitized_strands.append(strand)

        return sanitized_strands

    @staticmethod
    def strands_intersect(strand_1, strand_2):
        if strand_1['start_index'] > strand_2['start_index']:
            temp = strand_2
            strand_2 = strand_1
            strand_1 = temp

        start_1 = strand_1['start_index']
        start_2 = strand_2['start_index']
        domains_1 = strand_types[strand_1['strand_name']].domains
        domains_2 = strand_types[strand_2['strand_name']].domains
        diff = start_2 - start_1
        for i in range(diff, len(domains_1)):
            if domains_1[i] == domains_2[i-diff]:
                return True

        return False

    def svg_draw_contents(self, name=None, num_instructions=None):
        if self._dwg is None:
            name = name if name is not None else 'output'
            self._svg_vertical_offset = 55
            width = str(10 + (self._total_domains + 10) * self._svg_domain_length) + "mm"
            height = "100%" if num_instructions is None \
                else str(5 + (num_instructions + 1) * self._svg_vertical_offset) + "mm"
            self._dwg = svgwrite.Drawing(name + '.svg', size=(width, height))

        self._svg_draw_register_outline()
        self.svg_draw_strands(self.coverings, 1)

    def svg_increment_vertical_offset(self):
        self._svg_vertical_offset += 50

    def _svg_draw_register_outline(self):
        self._dwg.add(self._dwg.line(("5mm", str(self._svg_vertical_offset) + "mm"),
                                     (str(self._svg_left_offset + self._total_domains * self._svg_domain_length) + "mm",
                                      str(self._svg_vertical_offset) + "mm"),
                                     stroke=svgwrite.rgb(0, 0, 0)))

        domains = 0
        for cell in self.cells:
            cell_type = cell_types[cell]
            num_domains = len(cell_type.domains)
            self._dwg.add(self._dwg.line((str(self._svg_left_offset + domains) + "mm", str(self._svg_vertical_offset) + "mm"),
                                         (str(self._svg_left_offset + domains) + "mm",
                                          str(self._svg_vertical_offset - self._svg_cell_height) + "mm"),
                                         stroke=svgwrite.rgb(0, 0, 0)))

            for i in range(1, num_domains):
                self._dwg.add(self._dwg.line((str(self._svg_left_offset + domains + i * self._svg_domain_length) + "mm",
                                              str(self._svg_vertical_offset) + "mm"),
                                             (str(self._svg_left_offset + domains + i * self._svg_domain_length) + "mm",
                                              str(self._svg_vertical_offset - self._svg_domain_length) + "mm"),
                                             stroke=svgwrite.rgb(0, 0, 0)))

            domains += num_domains * self._svg_domain_length

        self._dwg.add(self._dwg.line((str(self._svg_left_offset + domains) + "mm", str(self._svg_vertical_offset) + "mm"),
                                     (str(self._svg_left_offset + domains) + "mm",
                                      str(self._svg_vertical_offset - self._svg_cell_height) + "mm"),
                                     stroke=svgwrite.rgb(0, 0, 0)))

    def svg_draw_strands(self, strand_set, layer):
        if strand_set is None:
            return

        upper_y = str(self._svg_vertical_offset - (layer + 1) * self._svg_domain_length) + "mm"
        y = str(self._svg_vertical_offset - layer * self._svg_domain_length) + "mm"
        previous_domains = 0
        for cell_name in self.cells:
            cell = cell_types[cell_name]
            for i in range(len(cell.domains)):
                left = self._svg_left_offset + (i + previous_domains) * self._svg_domain_length
                one_sixth = str(left + self._svg_domain_length // 6) + "mm"
                one_third = str(left + self._svg_domain_length // 3) + "mm"
                two_thirds = str(left + 55 * self._svg_domain_length // 100) + "mm"
                five_sixths = str(left + 5 * self._svg_domain_length // 6) + "mm"
                right = str(left + self._svg_domain_length) + "mm"
                short_right = str(left + 3 * self._svg_domain_length // 5) + "mm"
                short_left = str(left + self._svg_domain_length // 3) + "mm"
                left = str(left) + "mm"
                domain_coverings, orthogonal_coverings = self.get_coverings_at_domain_index(previous_domains + i,
                                                                                            include_orthogonal=True,
                                                                                            strand_set=strand_set)
                strand = strand_types[domain_coverings[0]['strand_name']] if len(domain_coverings) > 0 else \
                    strand_types[orthogonal_coverings[0]['strand_name']] if len(orthogonal_coverings) > 0 else None
                color = 'rgb(0, 0, 0)' if strand is None else convert_hex_to_rgb(strand.color)
                if len(orthogonal_coverings) >= 1:
                    point_right = True
                    for covering in orthogonal_coverings:
                        if covering['start_index'] == previous_domains + i:
                            point_right = False
                            break

                    if point_right:
                        self._dwg.add(self._dwg.line((left, y), (right, upper_y), stroke=color,
                                                     stroke_width="1mm"))
                        if strand is not None and not strand.is_complementary:
                            self._svg_draw_upper_right_arrow(int(right[:-2]), int(upper_y[:-2]), color)
                    else:
                        self._dwg.add(self._dwg.line((left, upper_y), (right, y), stroke=color,
                                                     stroke_width="1mm"))
                        if strand is not None and strand.is_complementary:
                            self._svg_draw_upper_left_arrow(int(left[:-2]), int(upper_y[:-2]), color)

                if len(domain_coverings) == 1:
                    index = previous_domains + i - domain_coverings[0]['start_index']
                    if strand.is_complementary:
                        if index == 0:
                            self._svg_draw_left_arrow(int(short_left[:-2]), int(y[:-2]), color)
                        else:
                            self._dwg.add(
                                self._dwg.line((one_sixth, y), (one_third, y),
                                               stroke=color, stroke_width="1mm"))
                            self._dwg.add(
                                self._dwg.line((two_thirds, y), (five_sixths, y),
                                               stroke=color, stroke_width="1mm"))
                    else:
                        if index == len(strand.domains) - 1:
                            self._dwg.add(
                                self._dwg.line((left, y), (short_right, y), stroke=color, stroke_width="1mm"))
                            self._svg_draw_right_arrow(int(short_right[:-2]), int(y[:-2]), color)
                        elif index == 0:
                            self._dwg.add(
                                self._dwg.line((short_left, y), (right, y), stroke=color, stroke_width="1mm"))
                        else:
                            self._dwg.add(
                                self._dwg.line((left, y), (right, y), stroke=color, stroke_width="1mm"))
                elif len(domain_coverings) > 1:
                    y1 = str(self._svg_vertical_offset - (layer - 1) * self._svg_domain_length) + "mm"
                    y2 = str(self._svg_vertical_offset - layer * self._svg_domain_length) + "mm"
                    next_strand = strand_types[domain_coverings[1]['strand_name']]
                    self._dwg.add(self._dwg.line((left, y2), (right, y1), stroke=color,
                                                 stroke_width="1mm"))
                    self._dwg.add(self._dwg.line((left, y1), (right, y2), stroke=convert_hex_to_rgb(next_strand.color),
                                                 stroke_width="1mm"))

            previous_domains += len(cell.domains)

        if len(strand_set) >= 1:
            last_covering = strand_set[-1]
            strand = strand_types[last_covering['strand_name']]
            color = convert_hex_to_rgb(strand.color)
            for i in range(previous_domains, last_covering['start_index'] + len(strand.domains)):
                left = self._svg_left_offset + i * self._svg_domain_length
                right = str(left + self._svg_domain_length) + "mm"
                left = str(left) + "mm"
                self._dwg.add(self._dwg.line((left, y), (right, upper_y), stroke=color,
                                             stroke_width="1mm"))
                if not strand.is_complementary:
                    self._svg_draw_upper_right_arrow(int(right[:-2]), int(upper_y[:-2]), color)

    def _svg_draw_right_arrow(self, tip_x, tip_y, color):
        right = tip_x * 3.7795
        left = (tip_x - 2 * self._svg_domain_length // 3) * 3.7795
        y = tip_y * 3.7795
        upper_y = (tip_y - self._svg_domain_length // 3) * 3.7795
        lower_y = (tip_y + self._svg_domain_length // 3) * 3.7795
        self._dwg.add(
            self._dwg.polygon(points=[(right, y), (left, upper_y), (left, lower_y)],
                              stroke=color, fill=color, stroke_width="1mm"))

    def _svg_draw_left_arrow(self, tip_x, tip_y, color):
        left = (tip_x + self._svg_domain_length // 2) * 3.7795
        right = (tip_x + 2 * self._svg_domain_length // 3 + self._svg_domain_length // 2) * 3.7795
        y = tip_y * 3.7795
        upper_y = (tip_y - self._svg_domain_length // 3) * 3.7795
        lower_y = (tip_y + self._svg_domain_length // 3) * 3.7795
        self._dwg.add(
            self._dwg.polygon(points=[(left, y), (right, lower_y), (right, upper_y)],
                              stroke=color, fill=color, stroke_width="1mm"))

    def _svg_draw_upper_right_arrow(self, tip_x, tip_y, color):
        x1 = tip_x * 3.7795
        y1 = tip_y * 3.7795
        x2 = (tip_x - 2 * self._svg_domain_length // 3) * 3.7795
        y2 = (tip_y + self._svg_domain_length // 3) * 3.7795
        x3 = (tip_x - self._svg_domain_length // 4) * 3.7795
        y3 = (tip_y + 2 * self._svg_domain_length // 3) * 3.7795
        self._dwg.add(
            self._dwg.polygon(points=[(x1, y1), (x2, y2), (x3, y3)],
                              stroke=color, fill=color, stroke_width="1mm"))

    def _svg_draw_upper_left_arrow(self, tip_x, tip_y, color):
        x1 = tip_x * 3.7795
        y1 = tip_y * 3.7795
        x2 = (tip_x + self._svg_domain_length // 4) * 3.7795
        y2 = (tip_y + 2 * self._svg_domain_length // 3) * 3.7795
        x3 = (tip_x + 2 * self._svg_domain_length // 3) * 3.7795
        y3 = (tip_y + self._svg_domain_length // 3) * 3.7795
        self._dwg.add(
            self._dwg.polygon(points=[(x1, y1), (x2, y2), (x3, y3)],
                              stroke=color, fill=color, stroke_width="1mm"))

    def save_svg(self):
        if self._dwg is not None:
            self._dwg.save()
            self._dwg = None

    @staticmethod
    def decode_json(cells, coverings):
        self = Register()
        self.cells = cells

        self.coverings = []
        for covering in coverings:
            self.coverings.append({'start_index': covering['start_index'],
                                   'strand_name': covering['strand_name']})

        for cell in self.cells:
            self._total_domains += len(cell_types[cell].domains)

        return self


def convert_hex_to_rgb(hex_rgb):
    red = int(hex_rgb[1:3], 16)
    green = int(hex_rgb[3:5], 16)
    blue = int(hex_rgb[5:7], 16)
    return 'rgb(%d, %d, %d)' % (red, green, blue)


def add_cell_type():
    name = input('Enter cell type name: ')
    domains = input('Enter domain names, separated by commas: ').split(',')
    if name not in cell_types.keys():
        cell_types[name] = Cell(domains)


def add_cells_to_register():
    global registers

    register_name = input('Enter register name: ')
    if register_name not in registers:
        registers[register_name] = Register()

    current_register = registers[register_name]

    cell_name = input('Enter the cell name: ')
    if cell_name not in cell_types.keys():
        print('No such cell exists')
        return

    copies = int(input('Enter the number of copies of each cell: '))
    coverings = input('Enter initial coverings for each cell copy: ')

    if coverings == '':
        coverings = []
    else:
        coverings = coverings.split(';')

    cell_size = len(cell_types[cell_name].domains)
    for _ in range(copies):
        current_register.add_cell(cell_name)
        for covering in coverings:
            data = covering.split(',')
            strand_type = data[0]
            offset = int(data[1])
            current_register.attempt_attachment(-cell_size + offset, strand_type)


def add_strand_type():
    global strand_types
    name = input('What is the strand name? ')
    domains = input('Enter domains separated by commas: ').split(',')
    is_complementary = input('Is the strand complementary to the top strand? Y or N? ')
    color = input('What is the color hex code? ')
    if color == '':
        color = '#000000'

    strand_types[name] = Strand(domains, True if is_complementary.lower() == 'y' else False, color)


def add_instruction():
    instruction = []
    num_instructions = int(input('Enter the number of instruction strands: '))
    for _ in range(num_instructions):
        strand_name = input('Enter the strand name: ')
        instruction.append(strand_name)

    instructions.append(instruction)


def run_simulation():
    register_copies = copy.deepcopy(registers)
    for register_key in register_copies.keys():
        print(register_key)
        register = register_copies[register_key]

        total_domains = 0
        for cell_name in register.cells:
            total_domains += len(cell_types[cell_name].domains)

        for inst_num in range(len(instructions)):
            if show_unused_instruction_strands:
                unattached_matches = []
            else:
                unattached_matches = None

            pre_instruction_register = copy.deepcopy(register)
            register.svg_draw_contents(register_key, len(instructions))
            inst = instructions[inst_num]
            new_strands = []
            for _ in range(len(inst)):  # Repeat in case some strands should take effect after another
                displacement_occurred = True
                while displacement_occurred:  # Repeat in case of cascades
                    new_attachments = []
                    for strand_name in inst:
                        for i in range(total_domains):
                            new_attachment = register.attempt_attachment(i, strand_name, unattached_matches)
                            if new_attachment is not None:
                                new_attachments.extend(new_attachment)

                    # do first round of displacements preserving the new strands
                    if len(new_attachments) > 0:
                        register.displace_strands(new_attachments)
                    else:
                        displacement_occurred = False

                    displaced_strands = register.displace_strands()
                    if displaced_strands == new_attachments:  # all new strands did not stably bind
                        displacement_occurred = False
                    else:
                        new_strands.extend([strand for strand in new_attachments if strand not in displaced_strands])

                    if unattached_matches is not None:
                        unattached_matches.extend(displaced_strands)

            print("Instruction", inst_num + 1)
            new_strands.sort(key=lambda x: x['start_index'])
            if unattached_matches is not None:
                unattached_matches = [strand for strand in unattached_matches if strand not in new_strands]
                unattached_matches = register.sanitize_unattached_strands(unattached_matches, new_strands)

            if (len(new_strands) == 0 and (
                    unattached_matches is None or len(unattached_matches) == 0)) and inst_num > 0:
                print('No changes\n')
            else:
                pre_instruction_register.print(new_strands, unattached_matches)
                register.svg_draw_strands(new_strands, 3)
                register.svg_draw_strands(unattached_matches, 6)
                print()

            register.svg_increment_vertical_offset()

            if step_by_step_simulation:
                input('Press Enter to continue')

        print("Final result")
        register.print()
        print()
        register.svg_draw_contents()
        register.save_svg()
        if step_by_step_simulation:
            input('Press Enter to continue')

        if keep_results:
            registers[register_key] = register


def save_data():
    filename = input("Enter filename: ")
    with open(filename, 'w') as file:
        json.dump({
            'cell_types': cell_types,
            'strand_types': strand_types,
            'registers': registers,
            'instructions': instructions
        }, file, indent=4, cls=ObjectEncoder)
        file.flush()


def toggle_step_by_step_simulation():
    global step_by_step_simulation
    step_by_step_simulation = not step_by_step_simulation


def toggle_keep_results():
    global keep_results
    keep_results = not keep_results


def toggle_show_unused_instruction_strands():
    global show_unused_instruction_strands
    show_unused_instruction_strands = not show_unused_instruction_strands


def exit_loop():
    global program_loop
    program_loop = False


def decode_json_dict(d, cls):
    decoded_dict = {}
    for key in d.keys():
        decoded_dict[key] = cls.decode_json(**d[key])

    return decoded_dict


def simd_simulator(args):
    if len(args) > 1:
        print('Loading saved data...')
        with open(args[1]) as file:
            data = json.load(file)
            global cell_types, registers, strand_types, instructions
            cell_types = decode_json_dict(data['cell_types'], Cell)
            strand_types = decode_json_dict(data['strand_types'], Strand)
            registers = decode_json_dict(data['registers'], Register)
            instructions = data['instructions']

    choice_dict = {'1': add_cell_type,
                   '2': add_cells_to_register,
                   '3': add_strand_type,
                   '4': add_instruction,
                   '5': run_simulation,
                   '6': save_data,
                   '7': toggle_step_by_step_simulation,
                   '8': toggle_keep_results,
                   '9': toggle_show_unused_instruction_strands,
                   '10': exit_loop}

    while program_loop:
        choice = input('''Enter one of the following options:
1 - Add cell type
2 - Add cells to register
3 - Add strand type
4 - Add instruction
5 - Run simulation
6 - Save data
7 - Turn step-by-step simulation ''' + ('off\n' if step_by_step_simulation else 'on\n') +
                       '''8 - ''' + (
                           'Don\'t keep results after simulation\n' if keep_results else 'Keep results after simulation\n') +
                       '''9 - ''' + ('Don\'t Show unused instruction strands\n' if show_unused_instruction_strands
                                     else 'Show unused instruction strands\n') +
                       '''10 - Exit

''')

        if choice in choice_dict:
            choice_dict[choice]()


if __name__ == '__main__':
    simd_simulator(sys.argv)
