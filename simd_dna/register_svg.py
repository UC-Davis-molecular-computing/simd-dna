import svgwrite
from functools import partial
from simd_dna.functions import *
from simd_dna.classes import *
from typing import Optional, List


class RegisterSVGDrawing:
    _domain_length = 5  # millimeters
    _normal_size_parameters = {
        'left_offset': 40,
        'cell_label_height_offset': 5,
        'cell_height': 40,
        'initial_vertical_offset': 55,
        'vertical_offset_increment': 50,
        'layer_offset': 5
    }
    _compressed_size_parameters = {
        'left_offset': 10,
        'cell_label_height_offset': 5,
        'cell_height': 15,
        'initial_vertical_offset': 30,
        'vertical_offset_increment': 25,
        'layer_offset': 3
    }

    def __init__(self, compress_svg_drawings: bool = False, draw_inert_instructions: bool = False) -> None:
        self._dwg = None
        self.compress_svg_drawings = compress_svg_drawings
        self.draw_inert_instructions = draw_inert_instructions
        self._current_size_parameters = self._normal_size_parameters
        self._vertical_offset = self._current_size_parameters['initial_vertical_offset']

    def initialize(self, register: Register,
                   name: Optional[str] = None,
                   num_instructions: Optional[int] = None) -> None:
        name = name if name is not None else 'output'
        if self.compress_svg_drawings:
            self._current_size_parameters = self._compressed_size_parameters
        else:
            self._current_size_parameters = self._normal_size_parameters
        self._vertical_offset = self._current_size_parameters['initial_vertical_offset']
        width = str(10 + (register.total_domains + 10) * self._domain_length) + "mm"
        height = "100%" if num_instructions is None \
            else str(self._current_size_parameters['initial_vertical_offset']
                     + self._vertical_offset
                     + (num_instructions * self._current_size_parameters['vertical_offset_increment'])) + "mm"
        self._dwg = svgwrite.Drawing(name + '.svg', size=(width, height))

    def draw_contents(self, register: Register, label: Optional[str] = None, draw_x: bool = False) -> None:
        self._draw_register_outline(register, label, draw_x)
        self.draw_strands(register, register.top_strands, 1)
        self._draw_cell_strand_labels(register)

    def increment_vertical_offset(self) -> None:
        self._vertical_offset += self._current_size_parameters['vertical_offset_increment']

    def _draw_register_outline(self, register: Register, label: str, draw_x: bool = False) -> None:
        if label is not None:
            self._dwg.add(
                self._dwg.text(label, x=[str(float(self._current_size_parameters['left_offset'] / 2)) + "mm"],
                               y=[str(float((self._vertical_offset - self._current_size_parameters[
                                   'cell_height'] / 2))) + "mm"],
                               fill=svgwrite.rgb(0, 0, 0),
                               style="text-anchor:middle;dominant-baseline:middle;font-size:20;"
                                     "font-family:sans-serif"))

        if draw_x:
            left = (self._current_size_parameters['left_offset'] - 4) / 2
            right = left + 4
            up = self._vertical_offset - 3 * self._current_size_parameters['cell_height'] / 10
            down = up + 4
            self._dwg.add(self._dwg.line(
                (str(left) + "mm", str(up) + "mm"), (str(right) + "mm", str(down) + "mm"),
                stroke=svgwrite.rgb(255, 0, 0),
                stroke_width="1mm"))
            self._dwg.add(self._dwg.line(
                (str(left) + "mm", str(down) + "mm"), (str(right) + "mm", str(up) + "mm"),
                stroke=svgwrite.rgb(255, 0, 0),
                stroke_width="1mm"))

        self._dwg.add(self._dwg.line(
            (str(self._current_size_parameters['left_offset']) + "mm", str(self._vertical_offset) + "mm"),
            (str(self._current_size_parameters[
                     'left_offset'] + register.total_domains * self._domain_length) + "mm",
             str(self._vertical_offset) + "mm"),
            stroke=svgwrite.rgb(0, 0, 0)))

        domains = 0
        for cell in register.cells:
            cell_type = register.cell_types[cell]
            num_domains = len(cell_type.domains)
            self._dwg.add(
                self._dwg.line((str(self._current_size_parameters['left_offset'] + domains) + "mm",
                                str(self._vertical_offset) + "mm"),
                               (str(self._current_size_parameters['left_offset'] + domains) + "mm",
                                str(self._vertical_offset - self._current_size_parameters[
                                    'cell_height']) + "mm"),
                               stroke=svgwrite.rgb(0, 0, 0)))

            for i in range(1, num_domains):
                self._dwg.add(self._dwg.line((str(
                    self._current_size_parameters['left_offset'] + domains + i * self._domain_length) + "mm",
                                              str(self._vertical_offset) + "mm"),
                                             (str(self._current_size_parameters[
                                                      'left_offset'] + domains + i * self._domain_length) + "mm",
                                              str(self._vertical_offset - self._domain_length) + "mm"),
                                             stroke=svgwrite.rgb(0, 0, 0)))

            domains += num_domains * self._domain_length

        self._dwg.add(
            self._dwg.line((str(self._current_size_parameters['left_offset'] + domains) + "mm",
                            str(self._vertical_offset) + "mm"),
                           (str(self._current_size_parameters['left_offset'] + domains) + "mm",
                            str(self._vertical_offset - self._current_size_parameters['cell_height']) + "mm"),
                           stroke=svgwrite.rgb(0, 0, 0)))

    def draw_strands(self, register: Register,
                     strand_set: List[TopStrand],
                     layer: int,
                     is_unattached_set: bool = False) -> None:
        if strand_set is None:
            return

        non_complementary_stroke_dasharray = "4,2" if is_unattached_set and self.compress_svg_drawings else "1,0"
        complementary_stroke_dasharray = "4,2" if is_unattached_set or not self.compress_svg_drawings else "1,0"

        layer_offset = self._current_size_parameters['layer_offset'] if layer != 1 else self._domain_length
        y = self._vertical_offset - layer * layer_offset
        upper_y = str(y - self._domain_length) + "mm"
        y = str(y) + "mm"
        diagonal_strand_offset = 0.1323
        y_diagonal_offset = str((float(y[:-2]) + diagonal_strand_offset)) + "mm"
        previous_domains = 0
        current_start = None
        current_strand = None
        crossover_start = None
        crossover_strand = None
        delayed_draw_arrow = None
        for cell_name in register.cells:
            cell = register.cell_types[cell_name]
            for i in range(len(cell.domains)):
                left = self._current_size_parameters['left_offset'] + (
                        i + previous_domains) * self._domain_length
                short_right = str(left + 3 * self._domain_length // 5) + "mm"
                short_left = str(left + self._domain_length // 3) + "mm"
                top_strands, orthogonal_top_strands = register.get_top_strands_at_domain_index(previous_domains + i,
                                                                                               include_orthogonal=True,
                                                                                               strand_set=strand_set)
                strand = register.strand_types[top_strands[0].strand_name] if len(top_strands) > 0 else \
                    register.strand_types[orthogonal_top_strands[0].strand_name] if len(orthogonal_top_strands) > 0 \
                    else None
                color = convert_hex_to_rgb('#000000') if strand is None \
                    else convert_hex_to_rgb(strand.color)
                if len(orthogonal_top_strands) >= 1 and crossover_start is None:
                    orthogonal_strand = orthogonal_top_strands[0]
                    orthogonal_color = convert_hex_to_rgb(register.strand_types[orthogonal_strand.strand_name].color)
                    point_right = orthogonal_strand.start_index != previous_domains + i
                    previous_left = self._current_size_parameters['left_offset'] + (
                            i - 1 + previous_domains) * self._domain_length
                    right = str(previous_left + 3 * self._domain_length // 5) + "mm"

                    if current_start is not None and orthogonal_strand.strand_name == current_strand:
                        if register.strand_types[orthogonal_strand.strand_name].is_complementary:
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
                        right_minus = str(float(right[:-2]) - 0.5 + self._domain_length
                                          + diagonal_strand_offset) + "mm"
                        previous_right_minus = str(float(right[:-2]) - 0.5 + diagonal_strand_offset) + "mm"
                        if register.strand_types[orthogonal_strand.strand_name].is_complementary:
                            self._dwg.add(
                                self._dwg.line((previous_right_minus, y_diagonal_offset), (right_minus, upper_y_offset),
                                               stroke=orthogonal_color,
                                               stroke_width="1mm", stroke_dasharray=complementary_stroke_dasharray))
                        else:
                            self._draw_upper_right_arrow(float(right_minus[:-2]) + diagonal_strand_offset,
                                                         float(upper_y[:-2]) + diagonal_strand_offset,
                                                         orthogonal_color)
                            self._dwg.add(
                                self._dwg.line((previous_right_minus, y_diagonal_offset), (right_minus, upper_y_offset),
                                               stroke=orthogonal_color,
                                               stroke_width="1mm",
                                               stroke_dasharray=non_complementary_stroke_dasharray))
                    else:
                        left_plus = self._current_size_parameters['left_offset'] + (
                                i + previous_domains) * self._domain_length \
                                    + 2 * diagonal_strand_offset
                        right_plus = str(left_plus + self._domain_length + self._domain_length // 3
                                         + diagonal_strand_offset) + "mm"
                        left_plus = str(left_plus) + "mm"
                        if register.strand_types[orthogonal_strand.strand_name].is_complementary:
                            self._draw_upper_left_arrow(float(left_plus[:-2]), float(upper_y_offset[:-2]),
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

                if len(top_strands) > 1 or crossover_start is not None:
                    if crossover_start is None:
                        crossover_start = short_left
                        crossover_strand = top_strands[1].strand_name
                        crossover_domain_count = 0

                    next_top_strands, orthogonal_top_strands = register.get_top_strands_at_domain_index(
                        previous_domains + i + 1,
                        include_orthogonal=True,
                        strand_set=strand_set)
                    crossover_domain_count += 1

                    if len(next_top_strands) + len(orthogonal_top_strands) <= 1:
                        first_color = register.strand_types[current_strand].color
                        second_color = register.strand_types[crossover_strand].color
                        left_diagonal_start = float(crossover_start[:-2]) - 0.5 + diagonal_strand_offset
                        right_diagonal_start = str(left_diagonal_start - diagonal_strand_offset) + "mm"
                        left_diagonal_start = str(left_diagonal_start) + "mm"
                        left_diagonal_end = float(short_right[:-2]) + 0.5
                        right_diagonal_end = str(left_diagonal_end - diagonal_strand_offset) + "mm"
                        left_diagonal_end = str(left_diagonal_end) + "mm"
                        top_y = str(float(y[:-2]) + diagonal_strand_offset - (
                                crossover_domain_count - 0.5) * self._domain_length) + "mm"
                        self._draw_horizontal_line(strand, current_start, crossover_start, y, first_color, False,
                                                   non_complementary_stroke_dasharray,
                                                   complementary_stroke_dasharray)
                        self._dwg.add(self._dwg.line((right_diagonal_start, top_y), (right_diagonal_end,
                                                                                     y_diagonal_offset),
                                                     stroke=second_color,
                                                     stroke_width="1mm"))
                        # Draw the upper right arrow after the right horizontal strand is drawn
                        delayed_draw_arrow = partial(self._draw_upper_right_arrow,
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
                elif len(top_strands) == 1:
                    top_strand = top_strands[0]
                    index = previous_domains + i - top_strand.start_index
                    if current_start is None:
                        current_start = short_left
                        current_strand = top_strand.strand_name
                        if strand.is_complementary and index == 0:
                            self._draw_left_arrow(float(current_start[:-2]), int(y[:-2]), color)
                    elif index == len(strand.domains) - 1:
                        self._draw_horizontal_line(strand, current_start, short_right, y, color, True,
                                                   non_complementary_stroke_dasharray,
                                                   complementary_stroke_dasharray)
                        current_start = None
                        current_strand = None
                        if delayed_draw_arrow is not None:
                            delayed_draw_arrow()
                            delayed_draw_arrow = None

            previous_domains += len(cell.domains)

        if len(strand_set) >= 1:
            last_top_strand = strand_set[-1]
            strand = register.strand_types[last_top_strand.strand_name]
            color = convert_hex_to_rgb(strand.color)
            previous_left = self._current_size_parameters['left_offset'] + (
                    previous_domains - 1) * self._domain_length
            right = str(previous_left + 3 * self._domain_length // 5) + "mm"
            right_minus = str(float(right[:-2]) - 0.5) + "mm"
            last_index = last_top_strand.start_index + len(strand.domains)
            last_right = str(float(right_minus[:-2]) + (last_index - previous_domains) * self._domain_length) + "mm"
            upper_y = str(float(y[:-2]) -
                          (last_index - previous_domains) * self._domain_length) + "mm"

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
                    self._draw_upper_right_arrow(float(last_right[:-2]), float(upper_y[:-2]), color)
                    self._dwg.add(self._dwg.line((right_minus, y), (last_right, upper_y), stroke=color,
                                                 stroke_width="1mm",
                                                 stroke_dasharray=non_complementary_stroke_dasharray))

    def _draw_horizontal_line(self, strand: Strand, current_start: str,
                              short_right: str, y: str,
                              color: str, draw_arrow_head: bool,
                              non_complementary_stroke_dasharray: str,
                              complementary_stroke_dasharray: str) -> None:
        if strand.is_complementary:
            self._dwg.add(
                self._dwg.line((current_start, y), (short_right, y), stroke=color,
                               stroke_width="1mm", stroke_dasharray=complementary_stroke_dasharray))
        else:
            if draw_arrow_head:
                self._draw_right_arrow(int(short_right[:-2]), int(y[:-2]), color)
            self._dwg.add(
                self._dwg.line((current_start, y), (short_right, y), stroke=color,
                               stroke_width="1mm", stroke_dasharray=non_complementary_stroke_dasharray))

    def _draw_cell_strand_labels(self, register: Register) -> None:
        previous_domains = 0
        for cell_name in register.cells:
            cell = register.cell_types[cell_name]
            labels = register.cell_types[cell_name].strand_labels
            is_match = [True for _ in labels]
            for i in range(len(labels)):
                label = labels[i]
                for strand in label['strands']:
                    if previous_domains + strand[0] < 0:
                        continue  # todo: handle negative indices

                    top_strands, orthogonal_top_strands \
                        = register.get_top_strands_at_domain_index(previous_domains + strand[0], True)
                    top_strands = list(filter(lambda d: d.start_index == previous_domains + strand[0],
                                              top_strands))
                    top_strands = [d.strand_name for d in top_strands]
                    orthogonal_top_strands = list(filter(lambda d: d.start_index == previous_domains + strand[0],
                                                         orthogonal_top_strands))
                    orthogonal_top_strands = [d.strand_name for d in orthogonal_top_strands]
                    if strand[1] not in top_strands and strand[1] not in orthogonal_top_strands:
                        is_match[i] = False
                        break

            labels = [label['label'] for (label, match) in zip(labels, is_match) if match]
            if len(labels) > 0:
                left = previous_domains * self._domain_length
                right = left + len(cell.domains) * self._domain_length
                x = ((left + right) / 2) + self._current_size_parameters['left_offset']
                x = str(x) + "mm"
                self._dwg.add(self._dwg.text(labels[0], x=[x],
                                             y=[str(float((self._vertical_offset +
                                                           self._current_size_parameters[
                                                               'cell_label_height_offset']))) + "mm"],
                                             fill=svgwrite.rgb(0, 0, 0),
                                             style="text-anchor:middle;dominant-baseline:middle;font-size:22;"
                                                   "font-family:sans-serif"))

            previous_domains += len(cell.domains)

    def _draw_right_arrow(self, tip_x: float, tip_y: float, color: str) -> None:
        right = tip_x * 3.7795
        left = (tip_x - self._domain_length / 3) * 3.7795
        y = tip_y * 3.7795
        upper_y = (tip_y - self._domain_length / 8) * 3.7795
        lower_y = (tip_y + self._domain_length / 8) * 3.7795
        self._dwg.add(
            self._dwg.polygon(points=[(right, y), (left, upper_y), (left, lower_y)],
                              stroke=color, fill=color, stroke_width="1mm"))

    def _draw_left_arrow(self, tip_x: float, tip_y: float, color: str) -> None:
        left = tip_x * 3.7795
        right = left + (self._domain_length / 3) * 3.7795
        y = tip_y * 3.7795
        upper_y = (tip_y - self._domain_length / 8) * 3.7795
        lower_y = (tip_y + self._domain_length / 8) * 3.7795
        self._dwg.add(
            self._dwg.polygon(points=[(left, y), (right, lower_y), (right, upper_y)],
                              stroke=color, fill=color, stroke_width="1mm"))

    def _draw_upper_right_arrow(self, tip_x: float, tip_y: float, color: str) -> None:
        x1 = tip_x * 3.7795
        y1 = tip_y * 3.7795
        x2 = (tip_x - self._domain_length // 2) * 3.7795
        y2 = (tip_y + self._domain_length // 3) * 3.7795
        x3 = (tip_x - self._domain_length // 4) * 3.7795
        y3 = (tip_y + self._domain_length // 2) * 3.7795
        self._dwg.add(
            self._dwg.polygon(points=[(x1, y1), (x2, y2), (x3, y3)],
                              stroke=color, fill=color, stroke_width="1mm"))

    def _draw_upper_left_arrow(self, tip_x: float, tip_y: float, color: str) -> None:
        x1 = tip_x * 3.7795
        y1 = tip_y * 3.7795
        x2 = (tip_x + self._domain_length // 4) * 3.7795
        y2 = (tip_y + self._domain_length // 2) * 3.7795
        x3 = (tip_x + self._domain_length // 2) * 3.7795
        y3 = (tip_y + self._domain_length // 3) * 3.7795
        self._dwg.add(
            self._dwg.polygon(points=[(x1, y1), (x2, y2), (x3, y3)],
                              stroke=color, fill=color, stroke_width="1mm"))

    def save_svg(self) -> None:
        if self._dwg is not None:
            self._dwg.save(pretty=True)
            self._dwg = None
