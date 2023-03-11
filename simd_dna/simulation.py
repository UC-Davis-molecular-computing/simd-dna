import copy
from simd_dna.classes import *


class Simulation:
    """This is an object that contains all settings and methods used in the SIMD||DNA simulation

    :param step_by_step_simulation: A boolean that causes the console program to print each instruction one at a time if
        set to True. Otherwise, all instructions are printed at once in each cycle.
    :param keep_results: If set to True, the results after applying all instructions will overwrite the original
        contents of the register. Otherwise, the results will be discarded. Discarding the results can be helpful when
        testing out new instructions.
    :param show_inert_instruction_strands: If set to True, inert instruction strands will be printed on the terminal,
        alongside the applicable instructions. Otherwise, only applicable instruction strands are shown.

    :ivar Mapping[str, Strand] strand_types: A dict of possible strand types in the simulation, mapping the strand name
        to the :class:`simd_dna.classes.Strand` instance.
    :ivar Mapping[str, Cell] cell_types: A dict of possible cell types in the simulation, mapping the cell name to the
        :class:`simd_dna.classes.Cell` instance.
    :ivar Mapping[str, Register] registers: A dict of registers present in the solution, mapping the register name to
        the :class:`simd_dna.classes.Register` instance.
    :ivar List[List[str]] instructions: A list of instructions to be applied. Each instruction is a list of strings,
        which are the names of strand types present in that instruction. The instructions are applied in the order of
        their indices, starting from 0.
    """

    def __init__(self, step_by_step_simulation: bool = False,
                 keep_results: bool = False,
                 show_inert_instruction_strands: bool = False) -> None:
        self.strand_types = {}
        self.cell_types = {}
        self.registers = {}
        self.instructions = []
        self.step_by_step_simulation = step_by_step_simulation
        self.keep_results = keep_results
        self.show_inert_instruction_strands = show_inert_instruction_strands

    def add_cell_type(self, name: str, domains: List[str]) -> None:
        """Adds a new :class:`simd_dna.classes.Cell` type to the simulation.

        :param name: The name of the new cell type.
        :param domains: A list of strings representing the names of the domains that comprise the new cell type.
        """
        if name in self.cell_types.keys():
            raise ValueError('Another cell type already has that name')
        else:
            self.cell_types[name] = Cell(domains)

    def add_cells_to_register(self, register_name: str,
                              cell_name: str,
                              top_strands: Optional[List[TopStrand]] = None,
                              copies: int = 1) -> None:
        """Adds new cells to a register, starting from the right side of the register (the bottom strand's 5' end.)

        :param register_name: The name of the :class:`simd_dna.classes.Register` to add cells to. If the register
            doesn't already exist in the simulation, a new one is created.
        :param cell_name: The name of the cell type to be added to the register.
        :param top_strands: A list of top strands present on the cells to be added. If None, the added cells will be
            bare.
        :param copies: An integer representing the number of cell copies to be added.
        """
        if cell_name not in self.cell_types:
            raise ValueError('No such cell exists')
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
        """Adds a new :class:`simd_dna.classes.Strand` type to the simulation.

        :param name: The name of the new strand type.
        :param domains: A list of strings representing the names of the domains that comprise the new strand type.
        :param is_complementary: A boolean determining if the strand is complementary to the top strand. If True, the
            strand is designed to bind to and remove top strands from the register, instead of binding to the register's
            bottom strand.
        :param color: A string hexadecimal color code representing the color of this strand when drawn in SVG.
        """
        self.strand_types[name] = Strand(domains, is_complementary, color)

    def add_instruction(self, instruction_strands: List[str]) -> None:
        """Adds a new instruction to the simulation.

        :param instruction_strands: A list of strings representing the strand names present in the instruction.
        """
        self.instructions.append(instruction_strands)

    def add_cell_strand_label(self, cell_name: str, coordinate_strand_pairs: List[List], string_label: str) -> None:
        """Adds a new strand label to a cell. See :class:`simd_dna.classes.Cell` for a detailed description of strand
        labels.

        :param cell_name: The name of the cell to add a strand label to.
        :param coordinate_strand_pairs: A list of pairs of coordinates and strand names, representing the locations of
            the leftmost domains of each strand.
        :param string_label: A string that will be printed underneath the cell in the SVG drawing if the top strands
            match the coordinates and strands indicated.
        """
        if cell_name not in self.cell_types:
            raise ValueError('No such cell exists')
        self.cell_types[cell_name].add_strand_label(coordinate_strand_pairs, string_label)

    def run_instruction(self, register_name: str,
                        inst_num: int) -> Tuple[Register, Register, List[TopStrand], Optional[List[TopStrand]]]:
        """Applies an instruction to a register.

        :param register_name: The name of the affected :class:`simd_dna.classes.Register`
        :param inst_num: The integer index of the applicable instruction.
        :return: A tuple describing the results of applying the instruction: the :class:`simd_dna.classes.Register`
            after applying the instruction, a copy of the :class:`simd_dna.classes.Register` before the instruction,
            the list of applicable instruction strands, and (optionally) the list of inert instruction strands. The
            function will not keep track of inert instruction strands if show_inert_instruction_strands is set to False.
        """
        if register_name not in self.registers.keys():
            raise ValueError('No such register exists')

        if inst_num < 0 or inst_num >= len(self.instructions):
            raise ValueError('Invalid instruction index')

        register = self.registers[register_name]
        # If keep_results is False, a copy of the register is made so that the original remains unaffected
        if not self.keep_results:
            register = copy.deepcopy(register)
        before_register = copy.deepcopy(register)
        total_domains = 0
        for cell_name in register.cells:
            total_domains += len(self.cell_types[cell_name].domains)

        if self.show_inert_instruction_strands:
            inert_matches = []
        else:
            inert_matches = None

        inst = self.instructions[inst_num]
        new_strands = []
        for _ in range(len(inst)):  # Repeat in case some strands should take effect after another
            displacement_occurred = True
            while displacement_occurred:  # Repeat in case of toehold exchanges/cascades
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
