from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Union, Tuple
import re
from json import JSONEncoder


class Strand:
    """This is a representation of a DNA strand in the SIMD||DNA model.

    :param domains: A list of strings representing the domains of the strand in left to right order
    :param is_complementary: A boolean that indicates whether the strand is complementary to the top strand of the
        register or not. A top complementary strand in the SIMD||DNA model has the 3' end on the left and the 5' end
        on the right.
    :param color: A hexadecimal string that represents the strand's color when drawn in an SVG file
    """

    def __init__(self, domains: List[str], is_complementary: bool, color: str = '#000000') -> None:
        self.domains = domains
        self.is_complementary = is_complementary
        if re.match('^#?[A-Fa-f0-9]{6}$', color) is not None:
            self.color = color
            if not self.color.startswith('#'):
                self.color = '#' + self.color
        else:
            self.color = '#000000'

    @staticmethod
    def decode_json(domains: List[str], is_complementary: bool, color: str = '#000000', **kwargs) -> Strand:
        """Decodes a JSON object and returns an instance of :class:`simd_dna.classes.Strand`

        :param domains: A list of strings corresponding to the domains field
        :param is_complementary: A boolean corresponding to the is_complementary field
        :param color:  A string corresponding to the color field
        :param kwargs: kwargs is placed to avoid throwing errors in the decode step if excess data is present in the
        JSON object. Any excess data is ignored.
        :return: A :class:`simd_dna.classes.Strand` object
        """
        self = Strand(domains, is_complementary, color)
        return self


@dataclass
class TopStrand:
    """This is a representation of a top strand currently attached to a register.\n
    Attributes:\n
    **start_index:** The index of the strand's leftmost domain on the register, starting from the register's
    leftmost domain at index 0\n
    **strand_name:** The name of the strand at start_index's location
    """
    start_index: int
    strand_name: str

    @staticmethod
    def decode_json(start_index: int, strand_name: str, **kwargs) -> TopStrand:
        """Decodes a JSON object and returns an instance of :class:simd_dna.classes.TopStrand

        :param start_index: An integer representing the register position of the top strand's leftmost domain
        :param strand_name: A string representing the name/type of the top strand
        :param kwargs: kwargs is placed to avoid throwing errors in the decode step if excess data is present in the
            JSON object. Any excess data is ignored.
        :return: A :class:`simd_dna.classes.TopStrand` object
        """
        return TopStrand(start_index, strand_name)


class Cell:
    """This is a representation of a cell in the SIMD||DNA model. A cell is a unit of data in the register, and is
    further subdivided into domains, which consist of a small number of nucleotides.

    :param domains: A list of strings representing the domains of the cell in left to right order
    :param strand_labels: A list of dictionaries that map strand patterns to string labels. A string label will be
        written underneath the cell in the SVG drawing if the strand pattern matches. Each dictionary has the following
        key-value pairs:\n
        **strands:** A 2D list, where the first index is an integer that represents the start index of a strand relative
        to the first domain of the cell, starting at 0. The second index is the string name of the strand type that
        should be present at that index for the pattern to match.\n
        **label:** A string that will be printed underneath the cell in the SVG if the strand pattern in 'strands'
        matches the cell's current contents.\n
        One example is the following dictionary:\n
        ``{'strands': [[0, 'Zero-first'], [3, 'Zero-second']], 'label': '0'}``\n
        This means that if a strand of type 'Zero-first' has its leftmost domain attached to domain 0 of the cell, and
        a strand of type 'Zero-second' has its leftmost domain attached to domain 3 of the cell, then the label string
        '0' will be written underneath that cell in the SVG drawing.
    """

    def __init__(self, domains: List[str], strand_labels: Optional[List[Dict]] = None) -> None:
        if strand_labels is None:
            strand_labels = []

        self.domains = domains
        self.strand_labels = strand_labels

    @staticmethod
    def decode_json(domains: List[str], strand_labels: Optional[List[Dict]] = None, **kwargs) -> Cell:
        """Decodes a JSON object and returns an instance of :class:`simd_dna.classes.Cell` .

        :param domains: A list of strings corresponding to the domains field
        :param strand_labels: A list of dictionaries corresponding to the strand_labels field
        :param kwargs: kwargs is placed to avoid throwing errors in the decode step if excess data is present in the
            JSON object. Any excess data is ignored.
        :return: A :class:`simd_dna.classes.Cell` object
        """
        if strand_labels is None:
            strand_labels = []

        self = Cell(domains, strand_labels)
        return self

    def add_strand_label(self, coordinate_strand_pairs: List, string_label: str) -> None:
        """Adds a new strand label to the cell type. See :class:`simd_dna.classes.Cell` for a detailed breakdown of the
        strand label's data structure.

        :param coordinate_strand_pairs: A list, where the first item is an integer that represents the start index
            of a strand relative to the first domain of the cell, starting at 0. The second item is the string name of
            the strand type that should be present at that index for the pattern to match.
        :param string_label: A string that will be printed underneath the cell in the SVG if the strand pattern in
            'strands' matches the cell's current contents.
        """
        coordinate_strand_pairs.sort(key=lambda x: x[0])
        label = {'strands': coordinate_strand_pairs, 'label': string_label}
        self.strand_labels.append(label)


class Register:
    """This is a representation of a register in the SIMD||DNA model. A register, in practice, is a long DNA strand
    attached to a magnetic bead. This long DNA strand is referred to as the \"bottom strand\", and information is
    stored and encoded through nick patterns in the attached \"top strands.\" Instruction strands are applied to a
    register, where the top strands are altered through DNA strand displacement. Waste products are washed away before
    the next instruction strands are applied.

    :param cell_types: A dictionary of :class:`simd_dna.classes.Cell` instances representing the possible cell types
        that can be part of this :class:`simd_dna.classes.Register` instance. The dictionary maps strings, which
        represent the cell name, to the actual :class:`simd_dna.classes.Cell` instance.
    :param strand_types: A dictionary of :class:`simd_dna.classes.Strand` instances representing the possible strand
        types that can be part of this :class:`simd_dna.classes.Register` instance. The dictionary maps strings, which
        represent the strand name, to the actual :class:`simd_dna.classes.Strand` instance.

    :ivar List[str] cells: A list of :class:`simd_dna.classes.Cell` type names that represent the cells that compose
        this :class:`simd_dna.classes.Register`
    :ivar List[simd_dna.classes.TopStrand] top_strands: A list of :class:`simd_dna.classes.TopStrand` s present on the
        register.
    :ivar int total_domains: The total number of domains in the register's bottom strand
    """

    def __init__(self, cell_types: Optional[Dict[str]] = None, strand_types: Optional[Dict[str]] = None) -> None:
        if cell_types is None:
            cell_types = []

        if strand_types is None:
            strand_types = []

        self.cell_types = cell_types
        self.strand_types = strand_types
        self.cells = []
        self.top_strands = []
        self.total_domains = 0

    def add_cell(self, cell_name: str) -> None:
        """Adds a cell to the right of the register.

        :param cell_name: A string representing the cell type to be added
        """
        if cell_name not in self.cell_types.keys():
            raise ValueError('Cell type does not exist')

        self.cells.append(cell_name)
        self.total_domains += len(self.cell_types[cell_name].domains)

    def get_cell_at_domain_index(self, domain_index: int) -> (Optional[str], int):
        """Returns the name of the cell type at the given domain index, starting from index 0, as well as the numerical
        offset relative to the beginning of the enclosing cell, starting at 0. For example, if a register has cells
        of type A, B, C in that order, where each cell type has 3 domains, then domains 0-2 will return A, 3-5 will
        return B, and 6-8 will return C. The domain at index 3 is the 0th domain in cell B, so an offset of 0 will be
        returned.

        :param domain_index: The integer index of the domain in the register.
        :return: A tuple containing a string that represents the cell type name (or None if the domain index exceeds
            the total domain length of the register), and the integer offset of that domain relative to the start index
            of its enclosing cell (0 if the domain index exceeds the total domain length.)
        """

        total_domains = 0
        for cell_name in self.cells:
            cell = self.cell_types[cell_name]
            if total_domains <= domain_index < total_domains + len(cell.domains):
                return cell, domain_index - total_domains
            total_domains += len(cell.domains)

        return None, 0

    def get_top_strands_at_domain_index(self, domain_index: int,
                                        include_orthogonal: bool = False,
                                        strand_set: Optional[List[TopStrand]] = None) \
            -> Union[Tuple[List[TopStrand], List[TopStrand]], List[TopStrand]]:
        """Returns the DNA top strand(s) present at a given domain index.

        :param domain_index: The integer index of the domain in the register.
        :param include_orthogonal: A boolean specifying whether a separate list of DNA strands that are orthogonal at
            this domain index should be returned.
        :param strand_set: A list of DNA top strands to be inspected. The register's top_strands instance variable will
            be used if None.
        :return: The list of DNA top strands attached to the provided domain index, and the list of DNA top strands
            with orthogonal domains hanging above the provided domain index if include_orthogonal is set to
            true
        """

        top_strands = []
        cell, offset = self.get_cell_at_domain_index(domain_index)
        if cell is None:
            return top_strands
        domain_label = cell.domains[offset]

        if include_orthogonal:
            orthogonal_top_strands = []

        if strand_set is None:
            strand_set = self.top_strands

        # For every DNA top strand in strand_set, check if domain_index is within its index range
        # from start_index to start_index + number of strand domains
        # If it's within range, check if the domains of the bottom and top strand match, and add the strand to
        # top_strands if so
        # If include_orthogonal is set to true, add the strand to orthogonal_top_strands if the domains don't match
        for top_strand in strand_set:
            start_index = top_strand.start_index
            strand = self.strand_types[top_strand.strand_name]
            if start_index <= domain_index < start_index + len(strand.domains):
                if domain_label == strand.domains[domain_index - start_index]:
                    top_strands.append(top_strand)
                elif include_orthogonal:
                    orthogonal_top_strands.append(top_strand)

        if include_orthogonal:
            return top_strands, orthogonal_top_strands
        else:
            return top_strands

    def attempt_attachment(self, domain_index: int,
                           strand_type: str,
                           unattached_matches: Optional[List[TopStrand]] = None) -> Optional[List[TopStrand]]:
        """Attempts to attach a copy of strand_type, with its leftmost domain placed on top of the specified
        domain_index if it's complementary to the bottom strand (strand_type.is_complementary is False.) If
        complementary to the top strand, detaches all top strands that bind to strand_type from the register.

        :param domain_index: The integer index of the register domain that the leftmost domain of strand_type
            will attempt to attach to (e.g. attach strand 'one_first' starting at domain index 56.)
            If strand_type is complementary to the top strand, this parameter is ignored.
        :param strand_type: The name of the :class:`simd_dna.classes.Strand` that will be attached. The name must be one
            of the keys in the Register's strand_types instance variable
        :param unattached_matches: A list of :class:`simd_dna.classes.TopStrand` instances that complement the domains
            underneath, but are inert because no open toeholds are available. If the current strand to be attached
            matches but is inert, it will be added to this list. If the caller isn't interested in getting the location
            of inert instruction strands, None can be provided.
        :return: A list containing :class:`simd_dna.classes.TopStrand` instances of new strands if attachment was
            successful, or None if no strands attached
        """
        if strand_type not in self.strand_types.keys():
            raise ValueError('Strand type does not exist')

        # Allow negative indexing
        if domain_index < 0:
            total_domains = 0
            for cell_name in self.cells:
                total_domains += len(self.cell_types[cell_name].domains)

            domain_index += total_domains
        strand = self.strand_types[strand_type]

        if strand.is_complementary:
            # If strand_type is complementary to the top strand, store all the strands that top_strand removes in
            # displaced_strands, and store TopStrand instances in displacing_strands, which note the domain indices
            # where strand_type successfully binds to the top strand

            displaced_strands = []
            displacing_strands = []
            for top_strand in self.top_strands:
                top_strand_domains = self.strand_types[top_strand.strand_name].domains
                is_match = True
                for i in range(len(top_strand_domains)):
                    if strand.domains[i] != top_strand_domains[i]:
                        is_match = False
                        break

                if is_match:
                    strand_start = top_strand.start_index
                    strand_end = strand_start + len(top_strand_domains)
                    for i in range(strand_start, strand_end):
                        top_strands_at_domain = self.get_top_strands_at_domain_index(i)
                        # must have at least one insecure domain
                        if len(top_strands_at_domain) > 1 or top_strand not in top_strands_at_domain:
                            displaced_strands.append(top_strand)
                            displacing_strands.append(TopStrand(strand_start, strand_type))
                            break

            if len(displaced_strands) > 0:
                self.top_strands = [x for x in self.top_strands if x not in displaced_strands]
                return displacing_strands
            elif unattached_matches is not None:
                # If the caller provides an unattached_matches list and the strand doesn't successfully bind at
                # domain_index, check if strand_type matches at least one domain. If so, it's considered inert and
                # is added to the list if not already present.

                matchings = 0
                for i in range(len(strand.domains)):
                    cell, offset = self.get_cell_at_domain_index(domain_index + i)
                    if cell is not None and cell.domains[offset] == strand.domains[i]:
                        matchings += 1

                if matchings >= 1:
                    new_strand = TopStrand(domain_index, strand_type)
                    if new_strand not in unattached_matches:
                        unattached_matches.append(new_strand)
                        unattached_matches.sort(key=lambda x: x.start_index)
        else:
            # If strand_type is complementary to the bottom strand, check if an open toehold is present when the strand
            # is placed according to domain_index. If so, check if the strand matches at least two domains along the
            # bottom strand, and return a list containing a new TopStrand instance if so. If no open toeholds are
            # present, but the caller provides an unattached_matches list and there are at least two domain matches,
            # the strand is considered inert and is added to the list if not already present.

            has_open_toehold = False
            for i in range(len(strand.domains)):
                top_strands = self.get_top_strands_at_domain_index(domain_index + i)
                if len(top_strands) == 0:
                    cell, offset = self.get_cell_at_domain_index(domain_index + i)
                    if cell is not None and cell.domains[offset] == strand.domains[i]:
                        has_open_toehold = True
                        break

            if has_open_toehold or unattached_matches is not None:
                # must have at least two matching domains to attach
                matchings = 0
                for i in range(len(strand.domains)):
                    cell, offset = self.get_cell_at_domain_index(domain_index + i)
                    if cell is not None and cell.domains[offset] == strand.domains[i]:
                        matchings += 1

                if matchings >= 2:
                    new_top_strand = TopStrand(domain_index, strand_type)
                    if has_open_toehold:
                        self.top_strands.append(new_top_strand)
                        self.top_strands.sort(key=lambda x: x.start_index)
                        return [new_top_strand]
                    else:
                        if new_top_strand not in unattached_matches:
                            unattached_matches.append(new_top_strand)
                            unattached_matches.sort(key=lambda x: x.start_index)
                        return None

        return None

    def displace_strands(self, excluded_strands: Optional[List[TopStrand]] = None) -> List[TopStrand]:
        """Simulates DNA strand displacement on the register, initiated by the new strands introduced by
        :func:`simd_dna.classes.Register.attempt_attachment`. By first attaching all possible new strands before
        displacing existing strands on the register, cooperative strand displacement can be simulated.

        :param excluded_strands: A list of :class:`simd_dna.classes.TopStrand` s that should be exempt from displacement
        :return: A list of :class:`simd_dna.classes.TopStrand` s that were displaced
        """
        if excluded_strands is None:
            excluded_strands = []

        displaced_strands = []
        filtered_top_strands = [x for x in self.top_strands if x not in excluded_strands]

        for top_strand in filtered_top_strands:
            strand_start = top_strand.start_index
            strand = self.strand_types[top_strand.strand_name]
            strand_end = strand_start + len(strand.domains)
            insecure_domains = 0
            # count the number of domains in this strand that are contested by other strands, or orthogonal
            for i in range(strand_start, strand_end):
                top_strands_at_domain = self.get_top_strands_at_domain_index(i)
                if len(top_strands_at_domain) > 1 or top_strand not in top_strands_at_domain:
                    insecure_domains += 1

            # if the strand is bound to the register by 0 or 1 uncontested domains, it is displaced
            if insecure_domains >= strand_end - strand_start - 1:
                displaced_strands.append(top_strand)

        if len(displaced_strands) > 0:
            self.top_strands = [x for x in self.top_strands if x not in displaced_strands]

        return displaced_strands

    def print(self, new_strands: Optional[List[TopStrand]] = None,
              unused_strands: Optional[List[TopStrand]] = None) -> None:
        """Prints the register's current contents on the terminal.

        :param new_strands: A list of :class:`simd_dna.classes.TopStrand` s that will displace the current strands,
            which will be printed one level above the current strands on the register.
        :param unused_strands: A list of :class:`simd_dna.classes.TopStrand` s that would've attached to the register but
            are inert, which will be printed two levels above the current strands on the register.
        """
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
                _, orthogonal_top_strands = \
                    self.get_top_strands_at_domain_index(previous_domains + i,
                                                         include_orthogonal=True)
                if len(orthogonal_top_strands) >= 1:
                    point_right = True
                    for top_strand in orthogonal_top_strands:
                        if top_strand.start_index == previous_domains + i:
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

        if len(self.top_strands) >= 1:
            last_top_strand = self.top_strands[-1]
            strand = self.strand_types[last_top_strand.strand_name]
            for _ in range(previous_domains, last_top_strand.start_index + len(strand.domains)):
                print('/', end='')

        print()

        previous_domains = 0
        print('|', end='')
        for cell_name in self.cells:
            cell = self.cell_types[cell_name]
            for i in range(len(cell.domains)):
                top_strands = self.get_top_strands_at_domain_index(previous_domains + i)
                if len(top_strands) == 0:
                    print('□', end='')
                elif len(top_strands) == 1:
                    top_strand = top_strands[0]
                    strand = self.strand_types[top_strand.strand_name]
                    index = previous_domains + i - top_strand.start_index
                    if index < len(strand.domains) - 1:
                        print('=', end='')
                    else:
                        print('>', end='')
                else:
                    print('x', end='')

            print('|', end='')
            previous_domains += len(cell.domains)

        print()

    def _print_floating_strands(self, strand_set: List[TopStrand]) -> None:
        # Helper function for print(), which prints strands that are a few layers above the register's bottom strand.
        # In the convention we use, these strands
        previous_domains = 0
        print('|', end='')
        for cell_name in self.cells:
            cell = self.cell_types[cell_name]
            for i in range(len(cell.domains)):
                domain_top_strands, orthogonal_top_strands = self.get_top_strands_at_domain_index(previous_domains + i,
                                                                                                  include_orthogonal
                                                                                                  =True,
                                                                                                  strand_set
                                                                                                  =strand_set)
                if len(orthogonal_top_strands) >= 1:
                    point_right = True
                    for top_strand in orthogonal_top_strands:
                        if top_strand.start_index == previous_domains + i:
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
            last_top_strand = strand_set[-1]
            strand = self.strand_types[last_top_strand.strand_name]
            for _ in range(previous_domains, last_top_strand.start_index + len(strand.domains)):
                print('/', end='')

        print()

        previous_domains = 0
        print('|', end='')
        for cell_name in self.cells:
            cell = self.cell_types[cell_name]
            for i in range(len(cell.domains)):
                top_strands = self.get_top_strands_at_domain_index(previous_domains + i, strand_set=strand_set)
                if len(top_strands) == 0:
                    print(' ', end='')
                elif len(top_strands) == 1:
                    top_strand = top_strands[0]
                    strand = self.strand_types[top_strand.strand_name]
                    index = previous_domains + i - top_strand.start_index
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

    def _print_empty_layer(self, blank_char: str = ' ') -> None:
        # Helper function for print(), printing only the vertical lines that separate each cell in the register
        previous_domains = 0
        print('|', end='')
        for cell_name in self.cells:
            cell = self.cell_types[cell_name]
            for i in range(len(cell.domains)):
                print(blank_char, end='')

            print('|', end='')
            previous_domains += len(cell.domains)

        print()

    def sanitize_inert_strands(self, inert_strands: List[TopStrand],
                               new_strands: List[TopStrand]) -> List[TopStrand]:
        """Removes inert strands that overlap with other inert strands or newly attached strands. Two strands are said
        to overlap if the have complementary domains in common on the register. This is done to reduce the cluttering
        in the visual representation of the register, in case the user chooses to print the register on the console or
        display its contents through an SVG representation.

        :param inert_strands: A list of inert :class:`simd_dna.classes.TopStrand` s after applying an instruction
        :param new_strands: A list of newly attached :class:`simd_dna.classes.TopStrand` s after applying an instruction
        :return: A list of inert :class:`simd_dna.classes.TopStrand` s after overlapping strands are removed
        """
        sanitized_strands = []
        for strand in reversed(inert_strands):
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
                sanitized_strands.insert(0, strand)

        return sanitized_strands

    def strands_intersect(self, strand_1: TopStrand, strand_2: TopStrand) -> bool:
        """Checks if two strands occupy the same domain location(s), where their domain segments at that location are
        complementary to the bottom strand. The two strands compete over that domain(s) if so.

        :param strand_1: The first :class:`simd_dna.classes.TopStrand` to compare
        :param strand_2: The second :class:`simd_dna.classes.TopStrand` to compare
        :return: True if the strands intersect, False otherwise
        """
        if strand_1.start_index > strand_2.start_index:
            temp = strand_2
            strand_2 = strand_1
            strand_1 = temp

        start_1 = strand_1.start_index
        start_2 = strand_2.start_index
        domains_1 = self.strand_types[strand_1.strand_name].domains
        domains_2 = self.strand_types[strand_2.strand_name].domains
        diff = start_2 - start_1
        for i in range(diff, len(domains_1)):
            if domains_1[i] == domains_2[i - diff]:
                return True

        return False

    @staticmethod
    def decode_json(cell_types: List[Cell], strand_types: List[Strand], cells, **kwargs) -> Register:
        """Decodes a JSON object and returns an instance of :class:`simd_dna.classes.Register`

        :param cell_types: A list of :class:`simd_dna.classes.Cell` types that can be present on this register
        :param strand_types: A list of :class:`simd_dna.classes.Strand` types that can be present on this register
        :param cells: A string list of cell names that compose this Register instance, going from left to right as the
            list index increases
        :param kwargs: kwargs is placed to avoid throwing errors in the decode step if excess data is present in the
            JSON object. Any excess data is ignored.
        :return: A :class:`simd_dna.classes.Register` object
        """

        self = Register(cell_types, strand_types)
        self.cells = cells

        # 'coverings' is old name of 'top_strands'
        if 'coverings' in kwargs.keys():
            top_strands = kwargs['coverings']
        elif 'top_strands' in kwargs.keys():
            top_strands = kwargs['top_strands']
        else:
            top_strands = []

        self.top_strands = []
        for top_strand in top_strands:
            self.top_strands.append(TopStrand.decode_json(**top_strand))

        self.top_strands.sort(key=lambda x: x.start_index)

        for cell in self.cells:
            self.total_domains += len(self.cell_types[cell].domains)

        return self


class ObjectEncoder(JSONEncoder):
    """
    A JSONEncoder subclass that allows the json.dump() function to get the dictionary encoding of an object
    """
    def default(self, o):
        return o.__dict__
