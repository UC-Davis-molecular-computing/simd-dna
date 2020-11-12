import copy
import itertools
import sys
import json
from json import JSONEncoder

program_loop = True
cell_types = {}
strands = {}
registers = {}
instructions = []


class ObjectEncoder(JSONEncoder):
    def default(self, o):
        return o.__dict__


class Strand:

    def __init__(self, name, domains, is_complementary):
        self.name = name
        self.domains = domains
        self.is_complementary = is_complementary

    @staticmethod
    def decode_json(name, domains, is_complementary):
        self = Strand(name, domains, is_complementary)
        return self


class Cell:

    def __init__(self, name, domains):
        self.name = name
        self.domains = domains

    @staticmethod
    def decode_json(name, domains):
        self = Cell(name, domains)
        return self


class Register:

    def __init__(self, name):
        self.name = name
        self.cells = []
        self.coverings = []

    def add_cells(self, cell_data):
        self.cells.append(cell_data)

    def get_cell(self, index):
        if index < 0:
            index = len(self.cells) + index

        if index < len(self.cells):
            return self.cells[index]
        else:
            return None

    def add_covering(self, cell_index, strand_type, offset):
        if cell_index < 0:
            cell_index = len(self.cells) + cell_index
        cell = self.get_cell(cell_index)
        strand = strands[strand_type]

        if strand.is_complementary:
            return

        # must have at least two matching domains to attach
        matches = 0
        for dom1, dom2 in itertools.zip_longest(strand.domains, cell.domains[offset:]):
            matches += 1 if dom1 == dom2 else 0

        # todo: handle displacement
        if matches >= 2:
            total_domains = 0
            for i in range(cell_index):
                total_domains += len(self.cells[i].domains)

            self.coverings.append({'start_index': total_domains + offset, 'strand': strand})
            self.coverings.sort(key=lambda x: x['start_index'])

    def cell_intersects_strand(self, cell_index, covering_index):
        if covering_index >= len(self.coverings):
            return False

        cell_start = 0
        for i in range(cell_index):
            cell_start += len(self.cells[i].domains)

        cell_end = cell_start + len(self.cells[cell_index].domains)
        covering_start = self.coverings[covering_index]['start_index']
        covering_end = covering_start + len(self.coverings[covering_index]['strand'].domains)

        return covering_end >= cell_start and covering_start <= cell_end

    def cell_matches_strand_domain(self, cell_index, covering_index, domain_index):
        domain1 = self.cells[cell_index].domains[domain_index]

        covering = self.coverings[covering_index]
        cell_start = 0
        for i in range(cell_index):
            cell_start += len(self.cells[i].domains)

        offset = cell_start - covering['start_index']
        if domain_index + offset < 0 or domain_index + offset >= len(covering['strand'].domains):
            return None

        strand_index = domain_index + offset
        domain2 = covering['strand'].domains[strand_index]
        is_last = strand_index + 1 == len(covering['strand'].domains)

        return (domain1 == domain2, is_last)

    def apply_instruction(self, instruction):
        pass

    @staticmethod
    def decode_json(name, cells, coverings):
        self = Register(name)
        self.cells = []
        for cell in cells:
            self.cells.append(Cell.decode_json(**cell))

        self.coverings = []
        for covering in coverings:
            self.coverings.append({'start_index': covering['start_index'],
                                   'strand': Strand.decode_json(**(covering['strand']))})
        return self


def add_cell_type():
    name = input('Enter cell type name: ')
    domains = input('Enter domain names, separated by commas: ').split(',')
    if name not in cell_types.keys():
        cell_types[name] = Cell(name, domains)


def add_cells_to_register():
    global registers

    register_name = input('Enter register name: ')
    if register_name not in registers:
        registers[register_name] = Register(register_name)

    current_register = registers[register_name]

    cell_name = input('Enter the cell name: ')
    copies = int(input('Enter the number of copies of each cell: '))

    for _ in range(copies):
        current_register.add_cells(cell_types[cell_name])
    coverings = input('Enter initial coverings for each cell copy: ')

    if coverings == "":
        return

    coverings = coverings.split(';')

    for i in range(-1, -len(coverings) - 1, -1):
        for covering in coverings:
            data = covering.split(',')
            strand_name = data[0]
            offset = int(data[1])
            current_register.add_covering(i, strand_name, offset)


def add_strand():
    global strands
    name = input('What is the strand name? ')
    domains = input('Enter domains separated by commas: ').split(',')
    is_complementary = input('Is the strand complementary to the top strand? Y or N? ')

    strands[name] = Strand(name, domains, True if is_complementary.lower() == 'y' else False)


def add_instruction():
    instruction = []
    instruction_strands = input('Enter the strand names, separated by commas: ').split(',')

    for strand in instruction_strands:
        instruction.append({'strand': strand})

    instructions.append(instruction)


def run_simulation():
    register_copies = copy.deepcopy(registers)
    for register in register_copies.values():
        print_register(register)
        print()
    input('Press any key to continue\n')


def print_register(register):
    cells = register.cells
    coverings = register.coverings

    print('|', end='')
    for cell in cells:
        for _ in cell.domains:
            print(' ', end='')

        print('|', end='')

    print()

    current_strand_index = 0
    previous_domains = 0
    print('|', end='')
    for cell_index in range(len(cells)):
        cell = cells[cell_index]
        intersecting_coverings = []
        for i in range(current_strand_index, len(coverings)):
            if register.cell_intersects_strand(cell_index, i):
                intersecting_coverings.append(i)
            else:
                break

        if len(intersecting_coverings) > 0:
            for domain_index in range(len(cell.domains)):
                matchings = list(map(lambda x: register.cell_matches_strand_domain(cell_index, x, domain_index),
                                     intersecting_coverings))
                mismatches = 0
                matches = 0
                for matching in matchings:
                    if matching is not None:
                        if matching[0]:
                            matches += 1
                        else:
                            mismatches += 1

                if mismatches >= 1 or matches >= 2:
                    print('↗', end='')
                else:
                    print(' ', end='')

            current_strand_index += len(intersecting_coverings) - 1
        else:
            for _ in cell.domains:
                print(' ', end='')

        print('|', end='')
        previous_domains += len(cell.domains)

    if len(coverings) > 1:
        last_covering = coverings[-1]
        for _ in range(previous_domains, last_covering['start_index'] + len(last_covering['strand'].domains)):
            print('↗', end='')

    print()

    current_strand_index = 0
    print('|', end='')
    for cell_index in range(len(cells)):
        cell = cells[cell_index]
        intersecting_coverings = []
        for i in range(current_strand_index, len(coverings)):
            if register.cell_intersects_strand(cell_index, i):
                intersecting_coverings.append(i)
            else:
                break

        if len(intersecting_coverings) > 0:
            for domain_index in range(len(cell.domains)):
                matchings = list(map(lambda x: register.cell_matches_strand_domain(cell_index, x, domain_index),
                                     intersecting_coverings))
                matches = 0
                for matching in matchings:
                    if matching is not None and matching[0]:
                        matches += 1

                if matches == 0:
                    print('□', end='')
                elif matches == 1:
                    for i in range(len(matchings)):
                        if matchings[i] is not None and matchings[i][0]:
                            is_last = matchings[i][1]
                            break

                    if is_last:
                        print('>', end='')
                    else:
                        print('=', end='')
                else:
                    print('⭜', end='')

            current_strand_index += len(intersecting_coverings) - 1
        else:
            for _ in cell.domains:
                print('□', end='')

        print('|', end='')

    print()


def save_data():
    filename = input("Enter filename: ")
    with open(filename, 'w') as file:
        json.dump({
            'cell_types': cell_types,
            'strands': strands,
            'registers': registers,
            'instructions': instructions
        }, file, indent=4, cls=ObjectEncoder)


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
            global cell_types, registers, strands, instructions
            cell_types = decode_json_dict(data['cell_types'], Cell)
            strands = decode_json_dict(data['strands'], Strand)
            registers = decode_json_dict(data['registers'], Register)
            instructions = data['instructions']

    choice_dict = {'1': add_cell_type,
                   '2': add_cells_to_register,
                   '3': add_strand,
                   '4': add_instruction,
                   '5': run_simulation,
                   '6': save_data,
                   '7': exit_loop}

    while program_loop:
        choice = input('''Enter one of the following options:
1 - Add cell type
2 - Add cells to register
3 - Add strand type
4 - Add instruction
5 - Run simulation
6 - Save data
7 - Exit

''')

        if choice in choice_dict:
            choice_dict[choice]()


if __name__ == '__main__':
    simd_simulator(sys.argv)
