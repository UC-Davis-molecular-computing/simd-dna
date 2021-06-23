from typing import Optional, List
from simd_dna.classes import *


class Simulation:

    def __init__(self, step_by_step_simulation=False, keep_results=False, show_unused_instruction_strands=False):
        self.strand_types = {}
        self.cell_types = {}
        self.registers = {}
        self.instructions = []
        self.step_by_step_simulation = step_by_step_simulation
        self.keep_results = keep_results
        self.show_unused_instruction_strands = show_unused_instruction_strands

    def add_cell_type(self, name, domains):
        if name not in self.cell_types.keys():
            self.cell_types[name] = Cell(domains)

    def add_cells_to_register(self, register_name, cell_name, coverings: Optional[List] = None, copies=1) -> None:
        if coverings is None:
            coverings = []
        if register_name not in self.registers:
            self.registers[register_name] = Register(self.cell_types, self.strand_types)

        current_register = self.registers[register_name]
        cell_size = len(self.cell_types[cell_name].domains)
        for _ in range(copies):
            current_register.add_cell(cell_name)
            for covering in coverings:
                strand_type = covering[0]
                offset = covering[1]
                current_register.attempt_attachment(-cell_size + offset, strand_type)

    def add_strand_type(self, name, domains, is_complementary=False, color='#000000'):
        self.strand_types[name] = Strand(domains, is_complementary, color)

    def add_instruction(self, instruction_strands):
        self.instructions.append(instruction_strands)

    def add_cell_strand_label(self, cell_name, coordinate_strand_pairs, string_label):
        self.cell_types[cell_name].add_strand_label(coordinate_strand_pairs, string_label)
