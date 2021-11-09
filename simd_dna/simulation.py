import copy
from typing import Optional, List, Tuple
from simd_dna.classes import *


class Simulation:

    def __init__(self, step_by_step_simulation: bool = False,
                 keep_results: bool = False,
                 show_unused_instruction_strands: bool = False) -> None:
        self.strand_types = {}
        self.cell_types = {}
        self.registers = {}
        self.instructions = []
        self.step_by_step_simulation = step_by_step_simulation
        self.keep_results = keep_results
        self.show_unused_instruction_strands = show_unused_instruction_strands

    def add_cell_type(self, name: str, domains: List[str]) -> None:
        if name not in self.cell_types.keys():
            self.cell_types[name] = Cell(domains)

    def add_cells_to_register(self, register_name: str,
                              cell_name: str,
                              top_strands: Optional[List[TopStrand]] = None,
                              copies: int = 1) -> None:
        if top_strands is None:
            top_strands = []
        if register_name not in self.registers:
            self.registers[register_name] = Register(self.cell_types, self.strand_types)

        current_register = self.registers[register_name]
        cell_size = len(self.cell_types[cell_name].domains)
        for _ in range(copies):
            current_register.add_cell(cell_name)
            for top_strand in top_strands:
                strand_type = top_strand[0]
                offset = top_strand[1]
                current_register.attempt_attachment(-cell_size + offset, strand_type)

    def add_strand_type(self, name: str, domains: List[str],
                        is_complementary: bool = False, color: str = '#000000') -> None:
        self.strand_types[name] = Strand(domains, is_complementary, color)

    def add_instruction(self, instruction_strands: List[str]) -> None:
        self.instructions.append(instruction_strands)

    def add_cell_strand_label(self, cell_name: str, coordinate_strand_pairs: List[List], string_label: str) -> None:
        self.cell_types[cell_name].add_strand_label(coordinate_strand_pairs, string_label)

    def run_instruction(self, register_name: str,
                        inst_num: int) -> Tuple[Register, Register, List[TopStrand], Optional[List[TopStrand]]]:
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
            inert_matches = []
        else:
            inert_matches = None

        inst = self.instructions[inst_num]
        new_strands = []
        for _ in range(len(inst)):  # Repeat in case some strands should take effect after another
            displacement_occurred = True
            while displacement_occurred:  # Repeat in case of cascades
                new_attachments = []
                for strand_name in inst:
                    for i in range(total_domains):
                        new_attachment = register.attempt_attachment(i, strand_name, inert_matches)
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

                if inert_matches is not None:
                    inert_matches.extend(displaced_strands)

        new_strands.sort(key=lambda x: x.start_index)
        if inert_matches is not None:
            inert_matches = [strand for strand in inert_matches if strand not in new_strands]
            inert_matches = register.sanitize_inert_strands(inert_matches, new_strands)

        return register, before_register, new_strands, inert_matches
