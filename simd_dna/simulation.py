import copy
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

    def run_instruction(self, register_name, inst_num):
        if register_name not in self.registers.keys():
            raise ValueError('No such register exists')

        if inst_num < 0 or inst_num >= len(self.instructions):
            raise ValueError('Invalid instruction index')

        register = self.registers[register_name]
        if not self.keep_results:
            register = copy.deepcopy(register)
        before_register = copy.deepcopy(register)
        total_domains = 0
        for cell_name in register.cells:
            total_domains += len(self.cell_types[cell_name].domains)

        if self.show_unused_instruction_strands:
            unattached_matches = []
        else:
            unattached_matches = None

        inst = self.instructions[inst_num]
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
                    new_strands.extend(
                        [strand for strand in new_attachments if strand not in displaced_strands])

                if unattached_matches is not None:
                    unattached_matches.extend(displaced_strands)

        new_strands.sort(key=lambda x: x['start_index'])
        if unattached_matches is not None:
            unattached_matches = [strand for strand in unattached_matches if strand not in new_strands]
            unattached_matches = register.sanitize_unattached_strands(unattached_matches, new_strands)

        return register, before_register, new_strands, unattached_matches
