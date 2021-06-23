import svgwrite
import re
from functools import partial
from simd_dna.functions import *
from json import JSONEncoder


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

    def add_strand_label(self, coordinate_strand_pairs, string_label):
        coordinate_strand_pairs.sort(key=lambda x: x[0])
        label = {'strands': coordinate_strand_pairs, 'label': string_label}
        self.strand_labels.append(label)


class Register:
    compress_svg_drawings = False
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

    def __init__(self, cell_types=None, strand_types=None):
        self.cell_types = cell_types
        self.strand_types = strand_types
        self.cells = []
        self.coverings = []
        self._dwg = None
        self._svg_current_size_parameters = self._svg_normal_size_parameters
        self._svg_vertical_offset = self._svg_current_size_parameters['initial_vertical_offset']
        self._total_domains = 0

    def add_cell(self, cell_name):
        self.cells.append(cell_name)
        self._total_domains += len(self.cell_types[cell_name].domains)

    def get_cell_at_domain_index(self, domain_index):
        total_domains = 0
        for cell_name in self.cells:
            cell = self.cell_types[cell_name]
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
            strand = self.strand_types[covering['strand_name']]
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
                total_domains += len(self.cell_types[cell_name].domains)

            domain_index += total_domains
        strand = self.strand_types[strand_type]

        if strand.is_complementary:
            displaced_strands = []
            displacing_strands = []
            for covering in self.coverings:
                top_strand = self.strand_types[covering['strand_name']]
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
            strand = self.strand_types[covering['strand_name']]
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
            cell = self.cell_types[cell_name]
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
            strand = self.strand_types[last_covering['strand_name']]
            for _ in range(previous_domains, last_covering['start_index'] + len(strand.domains)):
                print('/', end='')

        print()

        previous_domains = 0
        print('|', end='')
        for cell_name in self.cells:
            cell = self.cell_types[cell_name]
            for i in range(len(cell.domains)):
                coverings = self.get_coverings_at_domain_index(previous_domains + i)
                if len(coverings) == 0:
                    print('□', end='')
                elif len(coverings) == 1:
                    strand = self.strand_types[coverings[0]['strand_name']]
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
            cell = self.cell_types[cell_name]
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
            strand = self.strand_types[last_covering['strand_name']]
            for _ in range(previous_domains, last_covering['start_index'] + len(strand.domains)):
                print('/', end='')

        print()

        previous_domains = 0
        print('|', end='')
        for cell_name in self.cells:
            cell = self.cell_types[cell_name]
            for i in range(len(cell.domains)):
                coverings = self.get_coverings_at_domain_index(previous_domains + i, strand_set=strand_set)
                if len(coverings) == 0:
                    print(' ', end='')
                elif len(coverings) == 1:
                    strand = self.strand_types[coverings[0]['strand_name']]
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
            cell = self.cell_types[cell_name]
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

    def strands_intersect(self, strand_1, strand_2):
        if strand_1['start_index'] > strand_2['start_index']:
            temp = strand_2
            strand_2 = strand_1
            strand_1 = temp

        start_1 = strand_1['start_index']
        start_2 = strand_2['start_index']
        domains_1 = self.strand_types[strand_1['strand_name']].domains
        domains_2 = self.strand_types[strand_2['strand_name']].domains
        diff = start_2 - start_1
        for i in range(diff, len(domains_1)):
            if domains_1[i] == domains_2[i - diff]:
                return True

        return False

    def svg_initialize(self, name=None, num_instructions=None,):
        if self._dwg is None:
            name = name if name is not None else 'output'
            if Register.compress_svg_drawings:
                self._svg_current_size_parameters = self._svg_compressed_size_parameters
            self._svg_vertical_offset = self._svg_current_size_parameters['initial_vertical_offset']
            width = str(10 + (self._total_domains + 10) * self._svg_domain_length) + "mm"
            height = "100%" if num_instructions is None \
                else str(self._svg_current_size_parameters['initial_vertical_offset']
                         + self._svg_vertical_offset
                         + (num_instructions * self._svg_current_size_parameters['vertical_offset_increment'])) + "mm"
            self._dwg = svgwrite.Drawing(name + '.svg', size=(width, height))

    def svg_draw_contents(self, label=None):
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
            cell_type = self.cell_types[cell]
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

    def svg_draw_strands(self, strand_set, layer, is_unattached_set=False):
        if strand_set is None:
            return

        non_complementary_stroke_dasharray = "4,2" if is_unattached_set and Register.compress_svg_drawings else "1,0"
        complementary_stroke_dasharray = "4,2" if is_unattached_set or not Register.compress_svg_drawings else "1,0"

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
            cell = self.cell_types[cell_name]
            for i in range(len(cell.domains)):
                left = self._svg_current_size_parameters['left_offset'] + (
                        i + previous_domains) * self._svg_domain_length
                short_right = str(left + 3 * self._svg_domain_length // 5) + "mm"
                short_left = str(left + self._svg_domain_length // 3) + "mm"
                domain_coverings, orthogonal_coverings = self.get_coverings_at_domain_index(previous_domains + i,
                                                                                            include_orthogonal=True,
                                                                                            strand_set=strand_set)
                strand = self.strand_types[domain_coverings[0]['strand_name']] if len(domain_coverings) > 0 else \
                    self.strand_types[orthogonal_coverings[0]['strand_name']] if len(orthogonal_coverings) > 0 else None
                color = convert_hex_to_rgb('#000000') if strand is None \
                    else convert_hex_to_rgb(strand.color)
                if len(orthogonal_coverings) >= 1 and crossover_start is None:
                    orthogonal_strand = orthogonal_coverings[0]
                    orthogonal_color = convert_hex_to_rgb(self.strand_types[orthogonal_strand['strand_name']].color)
                    point_right = orthogonal_strand['start_index'] != previous_domains + i
                    previous_left = self._svg_current_size_parameters['left_offset'] + (
                            i - 1 + previous_domains) * self._svg_domain_length
                    right = str(previous_left + 3 * self._svg_domain_length // 5) + "mm"

                    if current_start is not None and orthogonal_strand['strand_name'] == current_strand:
                        if self.strand_types[orthogonal_strand['strand_name']].is_complementary:
                            self._dwg.add(
                                self._dwg.line((current_start, y), (right, y), stroke=orthogonal_color,
                                               stroke_width="1mm", stroke_dasharray=complementary_stroke_dasharray))
                        else:
                            self._dwg.add(
                                self._dwg.line((current_start, y), (right, y), stroke=orthogonal_color,
                                               stroke_width="1mm", stroke_dasharray=non_complementary_stroke_dasharray))
                        current_start = None
                        current_strand = None

                    upper_y_offset = str((float(upper_y[:-2]) + diagonal_strand_offset)) + "mm"
                    if point_right:
                        right_minus = str(float(right[:-2]) - 0.5 + self._svg_domain_length
                                          + diagonal_strand_offset) + "mm"
                        previous_right_minus = str(float(right[:-2]) - 0.5 + diagonal_strand_offset) + "mm"
                        if self.strand_types[orthogonal_strand['strand_name']].is_complementary:
                            self._dwg.add(
                                self._dwg.line((previous_right_minus, y_diagonal_offset), (right_minus, upper_y_offset),
                                               stroke=orthogonal_color,
                                               stroke_width="1mm", stroke_dasharray=complementary_stroke_dasharray))
                        else:
                            self._svg_draw_upper_right_arrow(float(right_minus[:-2]) + diagonal_strand_offset,
                                                             float(upper_y[:-2]) + diagonal_strand_offset,
                                                             orthogonal_color)
                            self._dwg.add(
                                self._dwg.line((previous_right_minus, y_diagonal_offset), (right_minus, upper_y_offset),
                                               stroke=orthogonal_color,
                                               stroke_width="1mm",
                                               stroke_dasharray=non_complementary_stroke_dasharray))
                    else:
                        left_plus = self._svg_current_size_parameters['left_offset'] + (
                                i + previous_domains) * self._svg_domain_length \
                                    + 2 * diagonal_strand_offset
                        right_plus = str(left_plus + self._svg_domain_length + self._svg_domain_length // 3 \
                                         + diagonal_strand_offset) + "mm"
                        left_plus = str(left_plus) + "mm"
                        if self.strand_types[orthogonal_strand['strand_name']].is_complementary:
                            self._svg_draw_upper_left_arrow(float(left_plus[:-2]), float(upper_y_offset[:-2]),
                                                            orthogonal_color)
                            self._dwg.add(self._dwg.line((left_plus, upper_y_offset), (right_plus, y_diagonal_offset),
                                                         stroke=orthogonal_color,
                                                         stroke_width="1mm",
                                                         stroke_dasharray=complementary_stroke_dasharray))
                        else:
                            self._dwg.add(self._dwg.line((left_plus, upper_y_offset), (right_plus, y_diagonal_offset),
                                                         stroke=orthogonal_color,
                                                         stroke_width="1mm",
                                                         stroke_dasharray=non_complementary_stroke_dasharray))

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
                        first_color = self.strand_types[current_strand].color
                        second_color = self.strand_types[crossover_strand].color
                        left_diagonal_start = float(crossover_start[:-2]) - 0.5 + diagonal_strand_offset
                        right_diagonal_start = str(left_diagonal_start - diagonal_strand_offset) + "mm"
                        left_diagonal_start = str(left_diagonal_start) + "mm"
                        left_diagonal_end = float(short_right[:-2]) + 0.5
                        right_diagonal_end = str(left_diagonal_end - diagonal_strand_offset) + "mm"
                        left_diagonal_end = str(left_diagonal_end) + "mm"
                        top_y = str(float(y[:-2]) + diagonal_strand_offset - (
                                crossover_domain_count - 0.5) * self._svg_domain_length) + "mm"
                        self._svg_draw_horizontal_line(strand, current_start, crossover_start, y, first_color, False,
                                                       non_complementary_stroke_dasharray,
                                                       complementary_stroke_dasharray)
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
                        self._svg_draw_horizontal_line(strand, current_start, short_right, y, color, True,
                                                       non_complementary_stroke_dasharray,
                                                       complementary_stroke_dasharray)
                        current_start = None
                        current_strand = None
                        if delayed_draw_arrow is not None:
                            delayed_draw_arrow()
                            delayed_draw_arrow = None

            previous_domains += len(cell.domains)

        if len(strand_set) >= 1:
            last_covering = strand_set[-1]
            strand = self.strand_types[last_covering['strand_name']]
            color = convert_hex_to_rgb(strand.color)
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
                                       stroke_width="1mm", stroke_dasharray=complementary_stroke_dasharray))
                else:
                    self._dwg.add(
                        self._dwg.line((current_start, y), (right, y), stroke=color,
                                       stroke_width="1mm", stroke_dasharray=non_complementary_stroke_dasharray))

            if last_index > previous_domains:
                if strand.is_complementary:
                    self._dwg.add(self._dwg.line((right_minus, y), (last_right, upper_y), stroke=color,
                                                 stroke_width="1mm", stroke_dasharray=complementary_stroke_dasharray))
                else:
                    self._svg_draw_upper_right_arrow(float(last_right[:-2]), float(upper_y[:-2]), color)
                    self._dwg.add(self._dwg.line((right_minus, y), (last_right, upper_y), stroke=color,
                                                 stroke_width="1mm",
                                                 stroke_dasharray=non_complementary_stroke_dasharray))

    def _svg_draw_horizontal_line(self, strand, current_start, short_right, y, color, draw_arrow_head,
                                  non_complementary_stroke_dasharray,
                                  complementary_stroke_dasharray):
        if strand.is_complementary:
            self._dwg.add(
                self._dwg.line((current_start, y), (short_right, y), stroke=color,
                               stroke_width="1mm", stroke_dasharray=complementary_stroke_dasharray))
        else:
            if draw_arrow_head:
                self._svg_draw_right_arrow(int(short_right[:-2]), int(y[:-2]), color)
            self._dwg.add(
                self._dwg.line((current_start, y), (short_right, y), stroke=color,
                               stroke_width="1mm", stroke_dasharray=non_complementary_stroke_dasharray))

    def _svg_draw_cell_strand_labels(self):
        previous_domains = 0
        for cell_name in self.cells:
            cell = self.cell_types[cell_name]
            labels = self.cell_types[cell_name].strand_labels
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
    def decode_json(cell_types, strand_types, cells, coverings, **kwargs):
        self = Register(cell_types, strand_types)
        self.cells = cells

        self.coverings = []
        for covering in coverings:
            self.coverings.append({'start_index': covering['start_index'],
                                   'strand_name': covering['strand_name']})

        for cell in self.cells:
            self._total_domains += len(self.cell_types[cell].domains)

        return self


class ObjectEncoder(JSONEncoder):
    def default(self, o):
        return o.__dict__
