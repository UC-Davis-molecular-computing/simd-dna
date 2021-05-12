import copy
import sys
import json
import svgwrite
import re
import ruamel.yaml as yaml
from json import JSONEncoder
from functools import partial
import colorsys

program_loop = True
step_by_step_simulation = False
keep_results = False
show_unused_instruction_strands = False
compress_svg_drawings = False
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
    def decode_json(domains, is_complementary, color='#000000', **kwargs):
        self = Strand(domains, is_complementary, color)
        return self


class Cell:

    def __init__(self, domains, strand_labels=[]):
        self.domains = domains
        self.strand_labels = strand_labels

    @staticmethod
    def decode_json(domains, strand_labels=[], **kwargs):
        self = Cell(domains, strand_labels)
        return self


class Register:
    _svg_domain_length = 5  # millimeters
    _svg_normal_size_parameters = {
        'left_offset': 40,
        'cell_label_height_offset': 5,
        'cell_height': 40,
        'initial_vertical_offset': 55,
        'vertical_offset_increment': 50,
        'layer_offset': 5
    }
    _svg_compressed_size_parameters = {
        'left_offset': 10,
        'cell_label_height_offset': 5,
        'cell_height': 15,
        'initial_vertical_offset': 30,
        'vertical_offset_increment': 25,
        'layer_offset': 3
    }

    def __init__(self):
        self.cells = []
        self.coverings = []
        self._dwg = None
        self._svg_current_size_parameters = self._svg_normal_size_parameters
        self._svg_vertical_offset = self._svg_current_size_parameters['initial_vertical_offset']
        self._total_domains = 0

    def add_cell(self, cell_name):
        self.cells.append(cell_name)
        self._total_domains += len(cell_types[cell_name].domains)

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

        sanitized_strands.sort(key=lambda x: x['start_index'])
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
            if domains_1[i] == domains_2[i - diff]:
                return True

        return False

    def svg_draw_contents(self, name=None, num_instructions=None, label=None):
        if self._dwg is None:
            name = name if name is not None else 'output'
            if compress_svg_drawings:
                self._svg_current_size_parameters = self._svg_compressed_size_parameters
            self._svg_vertical_offset = self._svg_current_size_parameters['initial_vertical_offset']
            width = str(10 + (self._total_domains + 10) * self._svg_domain_length) + "mm"
            height = "100%" if num_instructions is None \
                else str(self._svg_current_size_parameters['initial_vertical_offset']
                         + self._svg_vertical_offset
                         + (num_instructions * self._svg_current_size_parameters['vertical_offset_increment'])) + "mm"
            self._dwg = svgwrite.Drawing(name + '.svg', size=(width, height))

        self._svg_draw_register_outline(label)
        self.svg_draw_strands(self.coverings, 1)
        self._svg_draw_cell_strand_labels()

    def svg_increment_vertical_offset(self):
        self._svg_vertical_offset += self._svg_current_size_parameters['vertical_offset_increment']

    def _svg_draw_register_outline(self, label):
        if label is not None:
            self._dwg.add(
                self._dwg.text(label, x=[str(float(self._svg_current_size_parameters['left_offset'] / 2)) + "mm"],
                               y=[str(float((self._svg_vertical_offset - self._svg_current_size_parameters[
                                   'cell_height'] / 2))) + "mm"],
                               fill=svgwrite.rgb(0, 0, 0),
                               style="text-anchor:middle;dominant-baseline:middle;font-size:20;"
                                     "font-family:sans-serif"))

        self._dwg.add(self._dwg.line(
            (str(self._svg_current_size_parameters['left_offset']) + "mm", str(self._svg_vertical_offset) + "mm"),
            (str(self._svg_current_size_parameters[
                     'left_offset'] + self._total_domains * self._svg_domain_length) + "mm",
             str(self._svg_vertical_offset) + "mm"),
            stroke=svgwrite.rgb(0, 0, 0)))

        domains = 0
        for cell in self.cells:
            cell_type = cell_types[cell]
            num_domains = len(cell_type.domains)
            self._dwg.add(
                self._dwg.line((str(self._svg_current_size_parameters['left_offset'] + domains) + "mm",
                                str(self._svg_vertical_offset) + "mm"),
                               (str(self._svg_current_size_parameters['left_offset'] + domains) + "mm",
                                str(self._svg_vertical_offset - self._svg_current_size_parameters[
                                    'cell_height']) + "mm"),
                               stroke=svgwrite.rgb(0, 0, 0)))

            for i in range(1, num_domains):
                self._dwg.add(self._dwg.line((str(
                    self._svg_current_size_parameters['left_offset'] + domains + i * self._svg_domain_length) + "mm",
                                              str(self._svg_vertical_offset) + "mm"),
                                             (str(self._svg_current_size_parameters[
                                                      'left_offset'] + domains + i * self._svg_domain_length) + "mm",
                                              str(self._svg_vertical_offset - self._svg_domain_length) + "mm"),
                                             stroke=svgwrite.rgb(0, 0, 0)))

            domains += num_domains * self._svg_domain_length

        self._dwg.add(
            self._dwg.line((str(self._svg_current_size_parameters['left_offset'] + domains) + "mm",
                            str(self._svg_vertical_offset) + "mm"),
                           (str(self._svg_current_size_parameters['left_offset'] + domains) + "mm",
                            str(self._svg_vertical_offset - self._svg_current_size_parameters['cell_height']) + "mm"),
                           stroke=svgwrite.rgb(0, 0, 0)))

    def svg_draw_strands(self, strand_set, layer, use_lighter_colors=False):
        if strand_set is None:
            return

        layer_offset = self._svg_current_size_parameters['layer_offset'] if layer != 1 else self._svg_domain_length
        y = self._svg_vertical_offset - layer * layer_offset
        upper_y = str(y - self._svg_domain_length) + "mm"
        y = str(y) + "mm"
        diagonal_strand_offset = 0.1323
        y_diagonal_offset = str((float(y[:-2]) + diagonal_strand_offset)) + "mm"
        previous_domains = 0
        current_start = None
        current_strand = None
        crossover_start = None
        crossover_strand = None
        delayed_draw_arrow = None
        for cell_name in self.cells:
            cell = cell_types[cell_name]
            for i in range(len(cell.domains)):
                left = self._svg_current_size_parameters['left_offset'] + (
                            i + previous_domains) * self._svg_domain_length
                short_right = str(left + 3 * self._svg_domain_length // 5) + "mm"
                short_left = str(left + self._svg_domain_length // 3) + "mm"
                domain_coverings, orthogonal_coverings = self.get_coverings_at_domain_index(previous_domains + i,
                                                                                            include_orthogonal=True,
                                                                                            strand_set=strand_set)
                strand = strand_types[domain_coverings[0]['strand_name']] if len(domain_coverings) > 0 else \
                    strand_types[orthogonal_coverings[0]['strand_name']] if len(orthogonal_coverings) > 0 else None
                color = convert_hex_to_rgb('#000000', use_lighter_colors) if strand is None \
                    else convert_hex_to_rgb(strand.color, use_lighter_colors)
                if len(orthogonal_coverings) >= 1 and crossover_start is None:
                    orthogonal_strand = orthogonal_coverings[0]
                    orthogonal_color = convert_hex_to_rgb(strand_types[orthogonal_strand['strand_name']].color,
                                                          use_lighter_colors)
                    point_right = orthogonal_strand['start_index'] != previous_domains + i
                    previous_left = self._svg_current_size_parameters['left_offset'] + (
                                i - 1 + previous_domains) * self._svg_domain_length
                    right = str(previous_left + 3 * self._svg_domain_length // 5) + "mm"

                    if current_start is not None and orthogonal_strand['strand_name'] == current_strand:
                        if strand_types[orthogonal_strand['strand_name']].is_complementary:
                            self._dwg.add(
                                self._dwg.line((current_start, y), (right, y), stroke=orthogonal_color,
                                               stroke_width="1mm", stroke_dasharray="4,2"))
                        else:
                            self._dwg.add(
                                self._dwg.line((current_start, y), (right, y), stroke=orthogonal_color,
                                               stroke_width="1mm"))
                        current_start = None
                        current_strand = None

                    upper_y_offset = str((float(upper_y[:-2]) + diagonal_strand_offset)) + "mm"
                    if point_right:
                        right_minus = str(float(right[:-2]) - 0.5 + self._svg_domain_length
                                          + diagonal_strand_offset) + "mm"
                        previous_right_minus = str(float(right[:-2]) - 0.5 + diagonal_strand_offset) + "mm"
                        if strand_types[orthogonal_strand['strand_name']].is_complementary:
                            self._dwg.add(
                                self._dwg.line((previous_right_minus, y_diagonal_offset), (right_minus, upper_y_offset),
                                               stroke=orthogonal_color,
                                               stroke_width="1mm", stroke_dasharray="4,2"))
                        else:
                            self._svg_draw_upper_right_arrow(float(right_minus[:-2]) + diagonal_strand_offset,
                                                             float(upper_y[:-2]) + diagonal_strand_offset,
                                                             orthogonal_color)
                            self._dwg.add(
                                self._dwg.line((previous_right_minus, y_diagonal_offset), (right_minus, upper_y_offset),
                                               stroke=orthogonal_color,
                                               stroke_width="1mm"))
                    else:
                        left_plus = self._svg_current_size_parameters['left_offset'] + (
                                    i + previous_domains) * self._svg_domain_length \
                                    + 2 * diagonal_strand_offset
                        right_plus = str(left_plus + self._svg_domain_length + self._svg_domain_length // 3 \
                                         + diagonal_strand_offset) + "mm"
                        left_plus = str(left_plus) + "mm"
                        if strand_types[orthogonal_strand['strand_name']].is_complementary:
                            self._svg_draw_upper_left_arrow(float(left_plus[:-2]), float(upper_y_offset[:-2]),
                                                            orthogonal_color)
                            self._dwg.add(self._dwg.line((left_plus, upper_y_offset), (right_plus, y_diagonal_offset),
                                                         stroke=orthogonal_color,
                                                         stroke_width="1mm", stroke_dasharray="4,2"))
                        else:
                            self._dwg.add(self._dwg.line((left_plus, upper_y_offset), (right_plus, y_diagonal_offset),
                                                         stroke=orthogonal_color,
                                                         stroke_width="1mm"))

                if len(domain_coverings) > 1 or crossover_start is not None:
                    if crossover_start is None:
                        crossover_start = short_left
                        crossover_strand = domain_coverings[1]['strand_name']
                        crossover_domain_count = 0

                    next_coverings, orthogonal_coverings = self.get_coverings_at_domain_index(previous_domains + i + 1,
                                                                                              include_orthogonal=True,
                                                                                              strand_set=strand_set)
                    crossover_domain_count += 1

                    if len(next_coverings) + len(orthogonal_coverings) <= 1:
                        first_color = strand_types[current_strand].color
                        second_color = strand_types[crossover_strand].color
                        left_diagonal_start = float(crossover_start[:-2]) - 0.5 + diagonal_strand_offset
                        right_diagonal_start = str(left_diagonal_start - diagonal_strand_offset) + "mm"
                        left_diagonal_start = str(left_diagonal_start) + "mm"
                        left_diagonal_end = float(short_right[:-2]) + 0.5
                        right_diagonal_end = str(left_diagonal_end - diagonal_strand_offset) + "mm"
                        left_diagonal_end = str(left_diagonal_end) + "mm"
                        top_y = str(float(y[:-2]) + diagonal_strand_offset - (
                                crossover_domain_count - 0.5) * self._svg_domain_length) + "mm"
                        self._svg_draw_horizontal_line(strand, current_start, crossover_start, y, first_color, False)
                        self._dwg.add(self._dwg.line((right_diagonal_start, top_y), (right_diagonal_end,
                                                                                     y_diagonal_offset),
                                                     stroke=second_color,
                                                     stroke_width="1mm"))
                        # Draw the upper right arrow after the right horizontal strand is drawn
                        delayed_draw_arrow = partial(self._svg_draw_upper_right_arrow,
                                                     float(left_diagonal_end[:-2]), float(top_y[:-2]),
                                                     first_color)
                        self._dwg.add(
                            self._dwg.line((left_diagonal_start, y_diagonal_offset), (left_diagonal_end, top_y),
                                           stroke=first_color,
                                           stroke_width="1mm"))

                        current_start = short_right
                        current_strand = crossover_strand
                        crossover_start = None
                        crossover_strand = None
                elif len(domain_coverings) == 1:
                    index = previous_domains + i - domain_coverings[0]['start_index']
                    if current_start is None:
                        current_start = short_left
                        current_strand = domain_coverings[0]['strand_name']
                        if strand.is_complementary and index == 0:
                            self._svg_draw_left_arrow(float(current_start[:-2]), int(y[:-2]), color)
                    elif index == len(strand.domains) - 1:
                        self._svg_draw_horizontal_line(strand, current_start, short_right, y, color, True)
                        current_start = None
                        current_strand = None
                        if delayed_draw_arrow is not None:
                            delayed_draw_arrow()
                            delayed_draw_arrow = None

            previous_domains += len(cell.domains)

        if len(strand_set) >= 1:
            last_covering = strand_set[-1]
            strand = strand_types[last_covering['strand_name']]
            color = convert_hex_to_rgb(strand.color, use_lighter_colors)
            previous_left = self._svg_current_size_parameters['left_offset'] + (
                        previous_domains - 1) * self._svg_domain_length
            right = str(previous_left + 3 * self._svg_domain_length // 5) + "mm"
            right_minus = str(float(right[:-2]) - 0.5) + "mm"
            last_index = last_covering['start_index'] + len(strand.domains)
            last_right = str(float(right_minus[:-2]) + (last_index - previous_domains) * self._svg_domain_length) + "mm"
            upper_y = str(float(y[:-2]) -
                          (last_index - previous_domains) * self._svg_domain_length) + "mm"

            if current_start is not None:
                if strand.is_complementary:
                    self._dwg.add(
                        self._dwg.line((current_start, y), (right, y), stroke=color,
                                       stroke_width="1mm", stroke_dasharray="4,2"))
                else:
                    self._dwg.add(
                        self._dwg.line((current_start, y), (right, y), stroke=color,
                                       stroke_width="1mm"))

            if last_index > previous_domains:
                if strand.is_complementary:
                    self._dwg.add(self._dwg.line((right_minus, y), (last_right, upper_y), stroke=color,
                                                 stroke_width="1mm", stroke_dasharray="4,2"))
                else:
                    self._svg_draw_upper_right_arrow(float(last_right[:-2]), float(upper_y[:-2]), color)
                    self._dwg.add(self._dwg.line((right_minus, y), (last_right, upper_y), stroke=color,
                                                 stroke_width="1mm"))

    def _svg_draw_horizontal_line(self, strand, current_start, short_right, y, color, draw_arrow_head):
        if strand.is_complementary:
            self._dwg.add(
                self._dwg.line((current_start, y), (short_right, y), stroke=color,
                               stroke_width="1mm", stroke_dasharray="4,2"))
        else:
            if draw_arrow_head:
                self._svg_draw_right_arrow(int(short_right[:-2]), int(y[:-2]), color)
            self._dwg.add(
                self._dwg.line((current_start, y), (short_right, y), stroke=color,
                               stroke_width="1mm"))

    def _svg_draw_cell_strand_labels(self):
        global cell_types
        previous_domains = 0
        for cell_name in self.cells:
            cell = cell_types[cell_name]
            labels = cell_types[cell_name].strand_labels
            is_match = [True for _ in labels]
            for i in range(len(labels)):
                label = labels[i]
                for strand in label['strands']:
                    if previous_domains + strand[0] < 0:
                        continue  # todo: handle negative indices

                    domain_coverings, orthogonal_coverings \
                        = self.get_coverings_at_domain_index(previous_domains + strand[0], True)
                    domain_coverings = list(filter(lambda d: d['start_index'] == previous_domains + strand[0],
                                                   domain_coverings))
                    domain_coverings = [d['strand_name'] for d in domain_coverings]
                    orthogonal_coverings = list(filter(lambda d: d['start_index'] == previous_domains + strand[0],
                                                       orthogonal_coverings))
                    orthogonal_coverings = [d['strand_name'] for d in orthogonal_coverings]
                    if strand[1] not in domain_coverings and strand[1] not in orthogonal_coverings:
                        is_match[i] = False
                        break

            labels = [label['label'] for (label, match) in zip(labels, is_match) if match]
            if len(labels) > 0:
                left = previous_domains * self._svg_domain_length
                right = left + len(cell.domains) * self._svg_domain_length
                x = ((left + right) / 2) + self._svg_current_size_parameters['left_offset']
                x = str(x) + "mm"
                self._dwg.add(self._dwg.text(labels[0], x=[x],
                                             y=[str(float((self._svg_vertical_offset +
                                                           self._svg_current_size_parameters[
                                                               'cell_label_height_offset']))) + "mm"],
                                             fill=svgwrite.rgb(0, 0, 0),
                                             style="text-anchor:middle;dominant-baseline:middle;font-size:22;"
                                                   "font-family:sans-serif"))

            previous_domains += len(cell.domains)

    def _svg_draw_right_arrow(self, tip_x, tip_y, color):
        right = tip_x * 3.7795
        left = (tip_x - self._svg_domain_length / 3) * 3.7795
        y = tip_y * 3.7795
        upper_y = (tip_y - self._svg_domain_length / 8) * 3.7795
        lower_y = (tip_y + self._svg_domain_length / 8) * 3.7795
        self._dwg.add(
            self._dwg.polygon(points=[(right, y), (left, upper_y), (left, lower_y)],
                              stroke=color, fill=color, stroke_width="1mm"))

    def _svg_draw_left_arrow(self, tip_x, tip_y, color):
        left = tip_x * 3.7795
        right = left + (self._svg_domain_length / 3) * 3.7795
        y = tip_y * 3.7795
        upper_y = (tip_y - self._svg_domain_length / 8) * 3.7795
        lower_y = (tip_y + self._svg_domain_length / 8) * 3.7795
        self._dwg.add(
            self._dwg.polygon(points=[(left, y), (right, lower_y), (right, upper_y)],
                              stroke=color, fill=color, stroke_width="1mm"))

    def _svg_draw_upper_right_arrow(self, tip_x, tip_y, color):
        x1 = tip_x * 3.7795
        y1 = tip_y * 3.7795
        x2 = (tip_x - self._svg_domain_length // 2) * 3.7795
        y2 = (tip_y + self._svg_domain_length // 3) * 3.7795
        x3 = (tip_x - self._svg_domain_length // 4) * 3.7795
        y3 = (tip_y + self._svg_domain_length // 2) * 3.7795
        self._dwg.add(
            self._dwg.polygon(points=[(x1, y1), (x2, y2), (x3, y3)],
                              stroke=color, fill=color, stroke_width="1mm"))

    def _svg_draw_upper_left_arrow(self, tip_x, tip_y, color):
        x1 = tip_x * 3.7795
        y1 = tip_y * 3.7795
        x2 = (tip_x + self._svg_domain_length // 4) * 3.7795
        y2 = (tip_y + self._svg_domain_length // 2) * 3.7795
        x3 = (tip_x + self._svg_domain_length // 2) * 3.7795
        y3 = (tip_y + self._svg_domain_length // 3) * 3.7795
        self._dwg.add(
            self._dwg.polygon(points=[(x1, y1), (x2, y2), (x3, y3)],
                              stroke=color, fill=color, stroke_width="1mm"))

    def save_svg(self):
        if self._dwg is not None:
            self._dwg.save(pretty=True)
            self._dwg = None

    @staticmethod
    def decode_json(cells, coverings, **kwargs):
        self = Register()
        self.cells = cells

        self.coverings = []
        for covering in coverings:
            self.coverings.append({'start_index': covering['start_index'],
                                   'strand_name': covering['strand_name']})

        for cell in self.cells:
            self._total_domains += len(cell_types[cell].domains)

        return self


def convert_hex_to_rgb(hex_rgb, use_lighter_colors=False):
    red = int(hex_rgb[1:3], 16)
    green = int(hex_rgb[3:5], 16)
    blue = int(hex_rgb[5:7], 16)
    if use_lighter_colors:
        hue, lum, sat = colorsys.rgb_to_hls(red / 255, green / 255, blue / 255)
        lum = 0.8
        red, green, blue = tuple([int(color * 255) for color in colorsys.hls_to_rgb(hue, lum, sat)])
    return 'rgb(%d, %d, %d)' % (red, green, blue)


def convert_rgb_to_hex(r, g, b):
    red = int(r * 255)
    green = int(g * 255)
    blue = int(b * 255)
    return '#{:02x}{:02x}{:02x}'.format(red, green, blue)


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


def add_cell_strand_label():
    data = input('Enter the cell name, followed by the coordinate-strand name pairs, followed by the string label, \
all separated by commas: ').split(',')
    label = {'strands': [], 'label': data[-1]}
    for i in range(1, len(data) - 1, 2):
        strands = label['strands']
        index = int(data[i])
        if len(strands) == 0:
            strands.append([index, data[i + 1]])
        else:
            for j in range(len(strands)):
                if strands[j][0] > index:
                    strands.insert(j, [index, data[i + 1]])
                    break

                if j == len(strands) - 1:
                    strands.append([index, data[i + 1]])
    cell_types[data[0]].strand_labels.append(label)


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
            label = ("" if compress_svg_drawings else "Instruction ") + str(inst_num + 1)
            register.svg_draw_contents(register_key, len(instructions), label)
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
                register.svg_draw_strands(unattached_matches, 3 if compress_svg_drawings else 6, compress_svg_drawings)
                print()

            register.svg_increment_vertical_offset()

            if step_by_step_simulation:
                input('Press Enter to continue')

        print("Final result")
        register.print()
        print()
        label = "F" if compress_svg_drawings else "Final result"
        register.svg_draw_contents(register_key, label=label)
        register.save_svg()
        if step_by_step_simulation:
            input('Press Enter to continue')

        if keep_results:
            registers[register_key] = register


def save_data():
    filename = input("Enter filename: ")
    register_copies = copy.deepcopy(registers)
    for register in register_copies.values():
        del register._dwg
        del register._current_size_parameters
        del register._svg_vertical_offset
        del register._total_domains

    with open(filename, 'w') as file:
        json.dump({
            'cell_types': cell_types,
            'strand_types': strand_types,
            'registers': register_copies,
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


def toggle_compress_svg_drawings():
    global compress_svg_drawings
    compress_svg_drawings = not compress_svg_drawings


def convert_tm_to_simd():
    filename = input('Enter the YAML file name: ')
    with open(filename, 'r') as file:
        data = yaml.safe_load(file)
        blank = str(data['blank'])

        state_table_contains_unsupported_character = False
        transition_data = {}
        table = data['table']
        for state in table.keys():
            if state_table_contains_unsupported_character:
                break

            transition = table[state]
            if transition is not None:
                for key in transition.keys():
                    if 'write' in transition[key].keys() \
                            and str(transition[key]['write']) not in ['0', '1', blank]:
                        state_table_contains_unsupported_character = True
                        break

                    if type(key) is tuple:
                        if contains_outside_list(key, ['0', '1', blank]):
                            state_table_contains_unsupported_character = True
                            break
                        else:
                            for item in key:
                                add_simd_transition(state, str(item), copy.deepcopy(transition[key]), transition_data)
                    else:
                        if str(key) not in ['0', '1', blank]:
                            state_table_contains_unsupported_character = True
                            break
                        else:
                            add_simd_transition(state, str(key), transition[key], transition_data)

        if re.sub(rf'[{blank}01]', '', data['input']) != '' or state_table_contains_unsupported_character:
            print('Sorry, the TM conversion algorithm only supports the blank, 0, and 1 alphabet symbols.')
            return

        register_name = input('What name would you like to give the register? ')
        generate_tm_to_simd_data_from_transitions(transition_data, data, register_name)


def contains_outside_list(items, restricted_chars):
    for item in items:
        if str(item) not in restricted_chars:
            return True

    return False


def add_simd_transition(state, read, data, transition_data):
    key = (state, read)
    if 'write' not in data.keys():
        data['write'] = read
    else:
        data['write'] = str(data['write'])

    transition_data[key] = data


def generate_tm_to_simd_data_from_transitions(transition_data, tm_data, register_name):
    global cell_types, registers, strand_types, instructions
    cell_types = {}
    strand_types = {}
    registers = {}
    instructions = []

    # Create cell type
    tape_cell = Cell([])
    domain_template = '({},{})-{}'
    for key in transition_data.keys():
        for i in range(1, 3):
            tape_cell.domains.append(domain_template.format(*key, i))

    for i in range(1, 9):
        tape_cell.domains.append(str(i))

    cell_types['Tape cell'] = tape_cell

    create_basic_strand_types(transition_data, domain_template)
    create_tm_cell_labels(transition_data, tm_data['blank'])
    encode_register_data(register_name, transition_data, tm_data)
    generate_tm_instructions(transition_data, tm_data['blank'])


def create_basic_strand_types(transition_data, domain_template):
    global strand_types
    import seaborn as sns
    palette = sns.color_palette(None, len(transition_data) + 3)
    palette = [convert_rgb_to_hex(*x) for x in palette]

    strand_name_template = '({},{})_{}'
    strand_data = Strand([], False)
    for key, color in zip(transition_data.keys(), palette[:-3]):
        full_key = strand_name_template.format(*key, 'full')
        strand_types[full_key] = copy.deepcopy(strand_data)
        for i in range(1, 3):
            strand_types[full_key].domains.append(domain_template.format(*key, i))
            strand_types[full_key].color = color

    # Strand patterns
    # Blank - 123,456,78
    # Zero - 123,45,678
    # One - 12,345678

    blank_color = palette[-3]
    zero_color = palette[-2]
    one_color = palette[-1]
    strand_types['symbol_covered'] = copy.deepcopy(strand_data)
    for i in range(1, 9):
        strand_types['symbol_covered'].domains.append(str(i))

    strand_data.color = one_color
    strand_types['symbol_12'] = copy.deepcopy(strand_data)
    for i in range(1, 3):
        strand_types['symbol_12'].domains.append(str(i))

    strand_types['symbol_345678'] = copy.deepcopy(strand_data)
    for i in range(3, 9):
        strand_types['symbol_345678'].domains.append(str(i))

    strand_data.color = zero_color
    strand_types['symbol_123'] = copy.deepcopy(strand_data)
    for i in range(1, 4):
        strand_types['symbol_123'].domains.append(str(i))

    strand_types['symbol_45'] = copy.deepcopy(strand_data)
    for i in range(4, 6):
        strand_types['symbol_45'].domains.append(str(i))

    strand_types['symbol_678'] = copy.deepcopy(strand_data)
    for i in range(6, 9):
        strand_types['symbol_678'].domains.append(str(i))

    strand_data.color = blank_color
    strand_types['symbol_456'] = copy.deepcopy(strand_data)
    for i in range(4, 7):
        strand_types['symbol_456'].domains.append(str(i))

    strand_types['symbol_78'] = copy.deepcopy(strand_data)
    for i in range(7, 9):
        strand_types['symbol_78'].domains.append(str(i))


def create_tm_cell_labels(transition_data, blank_symbol):
    global cell_types
    label_template = {'strands': [], 'label': ''}
    current_index = 0
    for configuration in transition_data.keys():
        label_template['strands'].append([current_index, '({},{})_full'.format(*configuration)])
        current_index += 2

    blank_label = copy.deepcopy(label_template)
    blank_label['label'] = blank_symbol if blank_symbol != ' ' else '⎵'
    blank_label['strands'].append([current_index, 'symbol_123'])
    blank_label['strands'].append([current_index + 3, 'symbol_456'])
    blank_label['strands'].append([current_index + 6, 'symbol_78'])
    cell_types['Tape cell'].strand_labels.append(blank_label)

    zero_label = copy.deepcopy(label_template)
    zero_label['label'] = '0'
    zero_label['strands'].append([current_index, 'symbol_123'])
    zero_label['strands'].append([current_index + 3, 'symbol_45'])
    zero_label['strands'].append([current_index + 5, 'symbol_678'])
    cell_types['Tape cell'].strand_labels.append(zero_label)

    one_label = copy.deepcopy(label_template)
    one_label['label'] = '1'
    one_label['strands'].append([current_index, 'symbol_12'])
    one_label['strands'].append([current_index + 2, 'symbol_345678'])
    cell_types['Tape cell'].strand_labels.append(one_label)

    label_template['strands'].append([current_index, 'symbol_covered'])
    for i in range(len(transition_data)):
        data = list(transition_data.keys())[i]
        if data[1] == ' ':
            data = (data[0], '⎵')
        label_string = '({},{})'.format(*data)

        open_label = copy.deepcopy(label_template)
        open_label['strands'].pop(i)
        open_label['label'] = label_string
        cell_types['Tape cell'].strand_labels.append(open_label)

        plug_label = copy.deepcopy(label_template)
        plug_label['strands'][i][1] = plug_label['strands'][i][1].replace('full', 'plug')
        plug_label['label'] = label_string
        cell_types['Tape cell'].strand_labels.append(plug_label)

        final_label = copy.deepcopy(label_template)
        final_label['strands'][i][0] += -1
        final_label['strands'][i][1] = final_label['strands'][i][1].replace('full', 'final')
        final_label['label'] = label_string
        cell_types['Tape cell'].strand_labels.append(final_label)


def encode_register_data(register_name, transition_data, tm_data):
    global registers
    registers[register_name] = Register()
    if len(tm_data['input']) > 0:
        registers[register_name].add_cell('Tape cell')
        initial_symbol = tm_data['input'][0]
        initial_configuration = (tm_data['start state'], initial_symbol)
        current_index = 0
        has_valid_initial_transition = False
        for configuration in transition_data.keys():
            if initial_configuration == configuration:
                has_valid_initial_transition = True
            else:
                registers[register_name].coverings.append({
                    'start_index': current_index,
                    'strand_name': '({},{})_full'.format(*configuration)
                })

            current_index += 2
        if has_valid_initial_transition:
            registers[register_name].coverings.append({
                'start_index': current_index,
                'strand_name': 'symbol_covered'
            })
        elif initial_symbol == tm_data['blank']:
            insert_blank_symbol(registers[register_name].coverings, current_index)
        elif initial_symbol == '0':
            insert_zero_symbol(registers[register_name].coverings, current_index)
        else:
            insert_one_symbol(registers[register_name].coverings, current_index)
        current_index += 8
        for i in range(1, len(tm_data['input'])):
            symbol = tm_data['input'][i]
            registers[register_name].add_cell('Tape cell')
            for configuration in transition_data.keys():
                registers[register_name].coverings.append({
                    'start_index': current_index,
                    'strand_name': '({},{})_full'.format(*configuration)
                })
                current_index += 2
            if symbol == tm_data['blank']:
                insert_blank_symbol(registers[register_name].coverings, current_index)
            elif symbol == '0':
                insert_zero_symbol(registers[register_name].coverings, current_index)
            else:
                insert_one_symbol(registers[register_name].coverings, current_index)
            current_index += 8


def insert_blank_symbol(coverings, current_index):
    coverings.append({
        'start_index': current_index,
        'strand_name': 'symbol_123'
    })
    current_index += 3
    coverings.append({
        'start_index': current_index,
        'strand_name': 'symbol_456'
    })
    current_index += 3
    coverings.append({
        'start_index': current_index,
        'strand_name': 'symbol_78'
    })


def insert_zero_symbol(coverings, current_index):
    coverings.append({
        'start_index': current_index,
        'strand_name': 'symbol_123'
    })
    current_index += 3
    coverings.append({
        'start_index': current_index,
        'strand_name': 'symbol_45'
    })
    current_index += 2
    coverings.append({
        'start_index': current_index,
        'strand_name': 'symbol_678'
    })


def insert_one_symbol(coverings, current_index):
    coverings.append({
        'start_index': current_index,
        'strand_name': 'symbol_12'
    })
    current_index += 2
    coverings.append({
        'start_index': current_index,
        'strand_name': 'symbol_345678'
    })


def generate_tm_instructions(transition_data, blank_symbol):
    global strand_types, instructions
    configurations = list(transition_data.keys())
    for configuration in configurations:
        domains = []
        color = strand_types['({},{})_full'.format(*configuration)].color
        for i in range(1, 3):
            domains.append('({},{})-{}'.format(*configuration, i))
        domains.append('9')
        plug = Strand(domains, False, color)
        strand_types['({},{})_plug'.format(*configuration)] = plug
        unplug = copy.deepcopy(plug)
        unplug.is_complementary = True
        strand_types['({},{})_unplug'.format(*configuration)] = unplug

        domains = copy.deepcopy(domains)
        domains.pop()
        domains.insert(0, '10')
        final = Strand(domains, False, color)
        strand_types['({},{})_final'.format(*configuration)] = final
        restart = copy.deepcopy(final)
        restart.is_complementary = True
        strand_types['({},{})_restart'.format(*configuration)] = restart

    num_configurations = len(transition_data)
    for i in range(num_configurations):
        configuration = configurations[i]
        first_instruction = []
        if i == 0:
            for other_conf in configurations[1:]:
                first_instruction.append('({},{})_plug'.format(*other_conf))
        else:
            first_instruction.append('({},{})_unplug'.format(*configuration))

        if len(first_instruction) > 0:
            instructions.append(first_instruction)

        transition = transition_data[configuration]
        if 'L' in transition.keys():
            create_left_instruction(i, num_configurations, transition, transition_data, blank_symbol)
        else:
            create_right_instruction(i, num_configurations, transition, transition_data, blank_symbol)

    # Restart instruction
    instructions.append(['({},{})_restart'.format(*c) for c in configurations])


def create_left_instruction(index, num_configurations, transition, transition_data, blank_symbol):
    global strand_types, instructions
    next_state = transition['L']
    write_symbol = transition['write']
    configurations = list(transition_data.keys())

    # Displace all strands until left symbol area
    if index > 0:
        left_instruction = ['left_start']
        if 'left_start' not in strand_types.keys():
            first_domain = '({},{})-1'.format(*configurations[0])
            domains = ['7', '8', first_domain, '11']
            strand_types['left_start'] = Strand(domains, False)
            strand_types['left_start_strip'] = Strand(domains, True)
    else:
        left_instruction = ['left_start_extended']
        if 'left_start_extended' not in strand_types.keys():
            domains = ['7', '8']
            for i in range(1, 3):
                domains.append('({},{})-{}'.format(*configurations[0], i))
            domains.append('11')
            strand_types['left_start_extended'] = Strand(domains, False)
            strand_types['left_start_extended_strip'] = Strand(domains, True)
    for i in range(0, index - 1):
        current_configuration = configurations[i]
        next_configuration = configurations[i + 1]
        left_strand_name = '({},{})_left'.format(*current_configuration)
        if left_strand_name not in strand_types.keys():
            left_strand_domains = ['({},{})-2'.format(*current_configuration),
                                   '({},{})-1'.format(*next_configuration), '11']
            color = strand_types['({},{})_full'.format(*current_configuration)].color
            strand_types[left_strand_name] = Strand(left_strand_domains, False, color)
        left_instruction.append(left_strand_name)
    if index > 0:
        current_configuration = configurations[index - 1]
        next_configuration = configurations[index]
        left_extended_name = '({},{})_left_extended'.format(*current_configuration)
        if left_extended_name not in strand_types.keys():
            left_extended_domains = ['({},{})-2'.format(*current_configuration)]
            for i in range(1, 3):
                left_extended_domains.append('({},{})-{}'.format(*next_configuration, i))
            left_extended_domains.append('11')
            color = strand_types['({},{})_full'.format(*current_configuration)].color
            strand_types[left_extended_name] = Strand(left_extended_domains, False, color)
            left_extended_strip_name = '({},{})_left_extended_strip'.format(*current_configuration)
            strand_types[left_extended_strip_name] = Strand(left_extended_domains, True, color)
        left_instruction.append(left_extended_name)
    instructions.append(left_instruction)

    if 'one_left_start' not in strand_types.keys():
        first_domain = '({},{})-1'.format(*configurations[0])
        color = strand_types['symbol_345678'].color
        domains = ['3', '4', '5', '6', '7', '8', first_domain, '13']
        strand_types['one_left_start'] = Strand(domains, False, color)
        strand_types['one_left_strip'] = Strand(domains, True, color)

    # Add instructions for when left cell is zero
    next_index = configurations.index((next_state, '0')) \
        if (next_state, '0') in configurations else -1
    if next_index == -1:
        if 'temp_zero' not in strand_types.keys():
            color = strand_types['symbol_45'].color
            domains = ['4', '5', '6', '13']
            strand_types['temp_zero'] = Strand(domains, False, color)
            strand_types['temp_zero_strip'] = Strand(domains, True, color)
        instructions.append(['temp_zero'])
        instructions.append(['one_left_start'])
        instructions.append(['temp_zero_strip'])
        instructions.append(['symbol_45', 'symbol_678'])
    else:
        zero_instruction = []
        generate_left_cascade_instruction_strands(next_index, num_configurations, configurations, zero_instruction)
        if 'zero_symbol_displace' not in strand_types.keys():
            domains = ['2', '3', '4', '5', '6', '12']
            color = strand_types['symbol_45'].color
            strand_types['zero_symbol_displace'] = Strand(domains, False, color)
        zero_instruction.append('zero_symbol_displace')
        instructions.append(zero_instruction)
        instructions.append(['one_left_start'])  # in case left cell is 1
        generate_left_final_instruction_strands(next_index, num_configurations, configurations)

    # Add instructions for when left cell is one
    next_index = configurations.index((next_state, '1')) \
        if (next_state, '1') in configurations else -1
    instructions.append(['one_left_strip'])
    if next_index == -1:
        instructions.append(['symbol_345678'])
    else:
        if 'temp_one' not in strand_types.keys():
            color = strand_types['symbol_345678'].color
            first_domain = '({},{})-1'.format(*configurations[0])
            domains = ['4', '5', '6', '7', '8', first_domain, '13']
            strand_types['temp_one'] = Strand(domains, False, color)

        instructions.append(['temp_one'])

        one_instruction = []
        generate_left_cascade_instruction_strands(next_index, num_configurations, configurations, one_instruction)
        if 'one_symbol_displace' not in strand_types.keys():
            domains = ['2', '3', '12']
            color = strand_types['symbol_12'].color
            strand_types['one_symbol_displace'] = Strand(domains, False, color)
        one_instruction.append('one_symbol_displace')
        instructions.append(one_instruction)
        generate_left_final_instruction_strands(next_index, num_configurations, configurations)

    # Add instructions for when left cell is blank
    next_index = configurations.index((next_state, blank_symbol)) \
        if (next_state, blank_symbol) in configurations else -1
    if index > 0:
        instructions.append(['left_start_strip'])
    else:
        instructions.append(['left_start_extended_strip'])

    if next_index == -1:
        instructions.append(['symbol_78'])
    else:
        instructions.append(['symbol_right'])
        instructions.append(['symbol_right_strip'])
        blank_instruction = []
        generate_left_cascade_instruction_strands(next_index, num_configurations, configurations, blank_instruction)
        instructions.append(blank_instruction)
        generate_left_final_instruction_strands(next_index, num_configurations, configurations)

    # Replace current cell coverings
    if index > 0:
        instructions.append(['({},{})_left_extended_strip'.format(*configurations[index - 1])])
    else:
        instructions.append(['left_start_extended_strip'])

    replacement_instructions = []
    for i in range(0, index):
        replacement_instructions.append('({},{})_full'.format(*configurations[i]))
    generate_right_cascade_instruction_strands(index + 1, num_configurations, configurations, replacement_instructions)
    instructions.append(replacement_instructions)

    right_strip_instructions = []
    for i in range(index + 1, num_configurations):
        current_configuration = configurations[i]
        right_strip_instructions.append('({},{})_right_strip'.format(*current_configuration))
    right_strip_instructions.append('symbol_right_strip')
    instructions.append(right_strip_instructions)

    last_one_instructions = []
    for i in range(index, num_configurations):
        current_configuration = configurations[i]
        last_one_instructions.append('({},{})_full'.format(*current_configuration))
    if write_symbol == blank_symbol:
        last_one_instructions.extend(['symbol_123', 'symbol_456', 'symbol_78'])
    elif write_symbol == '0':
        last_one_instructions.extend(['symbol_123', 'symbol_45', 'symbol_678'])
    else:
        last_one_instructions.extend(['symbol_12', 'symbol_345678'])

    instructions.append(last_one_instructions)


def generate_left_cascade_instruction_strands(next_index, num_configurations, configurations, instruction_list):
    for i in range(next_index, num_configurations):
        current_configuration = configurations[i]
        next_configuration = configurations[i + 1] if i < num_configurations - 1 else -1
        left_strand_name = '({},{})_left'.format(*current_configuration)
        if left_strand_name not in strand_types.keys():
            left_strand_domains = ['({},{})-2'.format(*current_configuration)]
            if next_configuration == -1:
                left_strand_domains.append('1')
            else:
                left_strand_domains.append('({},{})-{}'.format(*next_configuration, 1))
            left_strand_domains.append('11')
            color = strand_types['({},{})_full'.format(*current_configuration)].color
            strand_types[left_strand_name] = Strand(left_strand_domains, False, color)
        instruction_list.append(left_strand_name)


def generate_left_final_instruction_strands(next_index, num_configurations, configurations):
    global instructions
    new_cell_instructions = ['({},{})_final'.format(*configurations[next_index])]
    for i in range(next_index + 1, num_configurations):
        new_cell_instructions.append('({},{})_full'.format(*configurations[i]))
    new_cell_instructions.append('symbol_covered')
    instructions.append(new_cell_instructions)


def create_right_instruction(index, num_configurations, transition, transition_data, blank_symbol):
    global strand_types, instructions
    next_state = transition['R']
    write_symbol = transition['write']
    configurations = list(transition_data.keys())

    # Displace all strands until right end of cell
    right_instruction = []

    generate_right_cascade_instruction_strands(index + 1, num_configurations, configurations, right_instruction)
    instructions.append(right_instruction)

    # Strip out strands from previous instruction
    right_strip_instructions = []
    for i in range(index + 1, num_configurations):
        current_configuration = configurations[i]
        right_strip_instructions.append('({},{})_right_strip'.format(*current_configuration))
    right_strip_instructions.append('symbol_right_strip')

    instructions.append(right_strip_instructions)

    # Replace current cell coverings
    replacement_instructions = []
    for current_configuration in configurations:
        replacement_instructions.append('({},{})_full'.format(*current_configuration))

    if write_symbol == blank_symbol:
        replacement_instructions.append('symbol_123')
        if 'symbol_45(13)' not in strand_types.keys():
            domains = ['4', '5', '13']
            strand_types['symbol_45(13)'] = Strand(domains, False, strand_types['symbol_456'].color)
            strand_types['symbol_45(13)_strip'] = Strand(domains, True, strand_types['symbol_456'].color)
            domains = ['6', '7', '13']
            strand_types['symbol_67(13)'] = Strand(domains, False, strand_types['symbol_78'].color)
            strand_types['symbol_67(13)_strip'] = Strand(domains, True, strand_types['symbol_78'].color)
        replacement_instructions.extend(['symbol_45(13)', 'symbol_67(13)'])
    elif write_symbol == '0':
        replacement_instructions.append('symbol_123')
        replacement_instructions.append('symbol_45')
        if 'symbol_67' not in strand_types.keys():
            strand_types['symbol_67'] = Strand(['6', '7'], False, strand_types['symbol_678'].color)
        replacement_instructions.append('symbol_67')
    else:
        replacement_instructions.append('symbol_12')
        if 'symbol_34567' not in strand_types.keys():
            strand_types['symbol_34567'] = Strand(['3', '4', '5', '6', '7'], False, strand_types['symbol_345678'].color)
        replacement_instructions.append('symbol_34567')

    instructions.append(replacement_instructions)

    # Displace right cell contents
    next_right_instructions = []
    if '({},{})_next_right'.format(*configurations[0]) not in strand_types.keys():
        for i in range(len(configurations)):
            domains = ['8' if i == 0 else '({},{})-2'.format(*configurations[i - 1])]
            current_configuration = configurations[i]
            domains.append('({},{})-1'.format(*current_configuration))
            domains.append('11')
            strand_name = '({},{})_next_right'.format(*current_configuration)
            color = strand_types['({},{})_full'.format(*current_configuration)].color
            strand_types[strand_name] = Strand(domains, False, color)
            strand_strip_name = '({},{})_next_right_strip'.format(*current_configuration)
            strand_types[strand_strip_name] = Strand(domains, True, color)

        last_configuration = configurations[num_configurations - 1]
        domains = ['({},{})-2'.format(*last_configuration), '1', '2', '11']
        strand_types['next_right_1'] = Strand(domains, False)
        strand_types['next_right_1_strip'] = Strand(domains, True)
        strand_types['next_right_2'] = Strand(['3', '4', '5', '11'], False)
        strand_types['next_right_2_strip'] = Strand(strand_types['next_right_2'].domains, True)

    for current_configuration in configurations:
        next_right_instructions.append('({},{})_next_right'.format(*current_configuration))
    next_right_instructions.append('next_right_1')
    next_right_instructions.append('next_right_2')

    instructions.append(next_right_instructions)

    # Add instructions for when right cell is blank
    next_index = configurations.index((next_state, blank_symbol)) \
        if (next_state, blank_symbol) in configurations else -1
    first_configuration = configurations[0]
    if 'right_blank' not in strand_types.keys():
        domains = ['4', '5', '6', '12']
        strand_types['right_blank'] = Strand(domains, False)
        strand_types['right_blank_strip'] = Strand(domains, True)
        domains = ['8']
        for i in range(1, 3):
            domains.append('({},{})-{}'.format(*first_configuration, i))
        domains.append('12')
        color = strand_types['({},{})_full'.format(*first_configuration)].color
        strand_types['8_right_plug'] = Strand(domains, False, color)
        strand_types['8_right_plug_strip'] = Strand(domains, True, color)
    blank_instructions = ['8_right_plug']
    for i in range(1, len(configurations)):
        current_configuration = configurations[i]
        strand_type = '({},{})_full' if i != next_index else '({},{})_final'
        blank_instructions.append(strand_type.format(*current_configuration))
    blank_instructions.append('symbol_123')
    blank_instructions.append('right_blank')
    instructions.append(blank_instructions)
    instructions.append(['right_blank_strip'])
    instructions.append(['symbol_456' if next_index == -1 else 'symbol_covered'])
    instructions.append(['8_right_plug_strip'])
    last_blank_instruction = []
    if write_symbol == blank_symbol:
        last_blank_instruction.extend(['symbol_456', 'symbol_78'])
    elif write_symbol == '0':
        last_blank_instruction.append('symbol_678')
    else:
        last_blank_instruction.append('symbol_345678')
    strand_type = '({},{})_final' if next_index == 0 else '({},{})_full'
    last_blank_instruction.append(strand_type.format(*first_configuration))
    instructions.append(last_blank_instruction)

    # Add instructions for when right cell is zero
    next_index = configurations.index((next_state, '0')) \
        if (next_state, '0') in configurations else -1
    instructions.append(['next_right_2_strip'])
    if next_index == -1:
        instructions.append(['symbol_123', 'symbol_45'])
    else:
        instructions.append(['symbol_covered'])
    zero_instructions = ['8_right_plug']
    for i in range(1, len(configurations)):
        current_configuration = configurations[i]
        strand_type = '({},{})_full' if i != next_index else '({},{})_final'
        zero_instructions.append(strand_type.format(*current_configuration))
    instructions.append(zero_instructions)
    instructions.append(['8_right_plug_strip'])
    last_zero_instruction = []
    if write_symbol == blank_symbol:
        last_zero_instruction.extend(['symbol_456', 'symbol_78'])
    elif write_symbol == '0':
        last_zero_instruction.append('symbol_678')
    else:
        last_zero_instruction.append('symbol_345678')
    strand_type = '({},{})_final' if next_index == 0 else '({},{})_full'
    last_zero_instruction.append(strand_type.format(*first_configuration))
    instructions.append(last_zero_instruction)

    # Add instructions for when right cell is one
    next_index = configurations.index((next_state, '1')) \
        if (next_state, '1') in configurations else -1
    one_instructions = []
    for current_configuration in configurations:
        one_instructions.append('({},{})_next_right_strip'.format(*current_configuration))
    one_instructions.append('next_right_1_strip')
    instructions.append(one_instructions)
    next_one_instruction = []
    if write_symbol == blank_symbol:
        next_one_instruction.extend(['symbol_456', 'symbol_78'])
    elif write_symbol == '0':
        next_one_instruction.append('symbol_678')
    else:
        next_one_instruction.append('symbol_345678')
    for i in range(len(configurations)):
        current_configuration = configurations[i]
        strand_type = '({},{})_full' if i != next_index else '({},{})_final'
        next_one_instruction.append(strand_type.format(*current_configuration))
    instructions.append(next_one_instruction)

    if next_index == -1:
        instructions.append(['symbol_12', 'symbol_345678'])
    else:
        instructions.append(['symbol_covered'])


def generate_right_cascade_instruction_strands(index, num_configurations, configurations, instruction_list):
    for i in range(index, num_configurations):
        previous_configuration = configurations[i - 1] if i > 0 else -1
        current_configuration = configurations[i]
        right_strand_name = '({},{})_right'.format(*current_configuration)
        if right_strand_name not in strand_types.keys():
            if previous_configuration == -1:
                right_strand_domains = ['8']
            else:
                right_strand_domains = ['({},{})-2'.format(*previous_configuration)]
            right_strand_domains.append('({},{})-1'.format(*current_configuration))
            right_strand_domains.append('11')
            color = strand_types['({},{})_full'.format(*current_configuration)].color
            strand_types[right_strand_name] = Strand(right_strand_domains, False, color)
            right_strip_name = '({},{})_right_strip'.format(*current_configuration)
            strand_types[right_strip_name] = Strand(right_strand_domains, True, color)
        instruction_list.append(right_strand_name)

    if 'symbol_right' not in strand_types.keys():
        last_configuration = configurations[num_configurations - 1]
        right_strand_domains = ['({},{})-2'.format(*last_configuration)]
        for i in range(1, 8):
            right_strand_domains.append(str(i))
        right_strand_domains.append('11')
        color = strand_types['symbol_covered'].color
        strand_types['symbol_right'] = Strand(right_strand_domains, False, color)
        strand_types['symbol_right_strip'] = Strand(right_strand_domains, True, color)

    instruction_list.append('symbol_right')


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
                   '5': add_cell_strand_label,
                   '6': run_simulation,
                   '7': save_data,
                   '8': toggle_step_by_step_simulation,
                   '9': toggle_keep_results,
                   '10': toggle_show_unused_instruction_strands,
                   '11': toggle_compress_svg_drawings,
                   '12': convert_tm_to_simd,
                   '13': exit_loop}

    while program_loop:
        choice = input('''Enter one of the following options:
1 - Add cell type
2 - Add cells to register
3 - Add strand type
4 - Add instruction
5 - Add cell-strand labels
6 - Run simulation
7 - Save data
8 - Turn step-by-step simulation ''' + ('off\n' if step_by_step_simulation else 'on\n') +
                       '''9 - ''' + (
                           'Don\'t keep results after simulation\n' if keep_results else 'Keep results after simulation\n') +
                       '''10 - ''' + ('Don\'t show unused instruction strands\n' if show_unused_instruction_strands
                                      else 'Show unused instruction strands\n') +
                       '''11 - ''' + ('Don\'t compress SVG drawings\n' if compress_svg_drawings
                                      else 'Compress SVG drawings\n') +
                       '''12 - Convert turingmachine.io Turing machine to SIMD register
13 - Exit

''')

        if choice in choice_dict:
            choice_dict[choice]()


if __name__ == '__main__':
    simd_simulator(sys.argv)
