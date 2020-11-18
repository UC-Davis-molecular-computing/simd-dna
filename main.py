import copy
import sys
import json
from json import JSONEncoder

program_loop = True
cell_types = {}
strand_types = {}
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

    def __init__(self):
        self.cells = []
        self.coverings = []

    def add_cell(self, cell_name):
        self.cells.append(cell_name)

    def get_cell_at_domain_index(self, domain_index):
        total_domains = 0
        for cell_name in self.cells:
            cell = cell_types[cell_name]
            if total_domains <= domain_index < total_domains + len(cell.domains):
                return cell, domain_index - total_domains
            total_domains += len(cell.domains)

        return None, 0

    def get_coverings_at_domain_index(self, domain_index):
        coverings = []
        cell, offset = self.get_cell_at_domain_index(domain_index)
        if cell is None:
            return coverings
        domain_label = cell.domains[offset]

        for covering in self.coverings:
            start_index = covering['start_index']
            strand = strand_types[covering['strand_name']]
            if start_index <= domain_index < start_index + len(strand.domains) \
                    and domain_label == strand.domains[domain_index - start_index]:
                coverings.append(covering)

        return coverings

    def attempt_displacement(self, domain_index, strand_type):
        if domain_index < 0:
            total_domains = 0
            for cell_name in self.cells:
                total_domains += len(cell_types[cell_name].domains)

            domain_index += total_domains
        strand = strand_types[strand_type]

        if strand.is_complementary:
            pass
        else:
            has_open_toehold = False
            for i in range(len(strand.domains)):
                coverings = self.get_coverings_at_domain_index(domain_index + i)
                if len(coverings) == 0:
                    cell, offset = self.get_cell_at_domain_index(domain_index + i)
                    if cell is not None and cell.domains[offset] == strand.domains[i]:
                        has_open_toehold = True

            if has_open_toehold:
                # must have at least two matching domains to attach
                matchings = 0
                for i in range(len(strand.domains)):
                    cell, offset = self.get_cell_at_domain_index(domain_index + i)
                    if cell is not None and cell.domains[offset] == strand.domains[i]:
                        matchings += 1

                if matchings >= 2:
                    new_covering = {'start_index': domain_index, 'strand_name': strand_type}
                    self.displace_strands(new_covering)

    def displace_strands(self, new_covering):
        self.coverings.append(new_covering)
        new_strand_start = self.coverings[-1]['start_index']
        new_strand = strand_types[self.coverings[-1]['strand_name']]
        new_strand_end = new_strand_start + len(new_strand.domains)
        displaced_strands = []
        for covering in self.coverings[:-1]:
            strand_start = covering['start_index']
            strand = strand_types[covering['strand_name']]
            strand_end = strand_start + len(strand.domains)
            if new_strand_end >= strand_start and new_strand_start < strand_end:
                insecure_domains = 0
                for i in range(strand_start, strand_end):
                    coverings_at_domain = self.get_coverings_at_domain_index(i)
                    if len(coverings_at_domain) > 1 or covering not in coverings_at_domain:
                        insecure_domains += 1

                if insecure_domains >= strand_end - strand_start - 1:
                    displaced_strands.append(covering)

        if len(displaced_strands) > 0:
            self.coverings = [x for x in self.coverings if x not in displaced_strands]

        # check if the new strand still stably binds
        insecure_domains = 0
        for i in range(new_strand_start, new_strand_end):
            coverings_at_domain = self.get_coverings_at_domain_index(i)
            if len(coverings_at_domain) > 1 or self.coverings[-1] not in coverings_at_domain:
                insecure_domains += 1

        if insecure_domains >= new_strand_end - new_strand_start - 1:
            self.coverings.pop(-1)
            displacement_succeeded = False
        else:
            displacement_succeeded = True

        self.coverings.sort(key=lambda x: x['start_index'])
        return displacement_succeeded

    @staticmethod
    def decode_json(cells, coverings):
        self = Register()
        self.cells = cells

        self.coverings = []
        for covering in coverings:
            self.coverings.append({'start_index': covering['start_index'],
                                   'strand_name': covering['strand_name']})
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
            current_register.attempt_displacement(-cell_size + offset, strand_type)


def add_strand_type():
    global strand_types
    name = input('What is the strand name? ')
    domains = input('Enter domains separated by commas: ').split(',')
    is_complementary = input('Is the strand complementary to the top strand? Y or N? ')

    strand_types[name] = Strand(name, domains, True if is_complementary.lower() == 'y' else False)


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
        print_register(register)
        print()
        input('Press any key to continue')

        total_domains = 0
        for cell_name in register.cells:
            total_domains += len(cell_types[cell_name].domains)

        for inst in instructions:
            for _ in range(len(inst)):  # Repeat in case some strands should take effect after another
                for strand_name in inst:
                    displacement_succeeded = True
                    while displacement_succeeded:
                        for i in range(total_domains):
                            success = register.attempt_displacement(i, strand_name)
                            if success:
                                break

                            if i == total_domains - 1:
                                displacement_succeeded = False

            print_register(register)
            print()
            input('Press any key to continue')


def print_register(register):
    cells = register.cells
    coverings = register.coverings

    print('|', end='')
    for cell_name in cells:
        cell = cell_types[cell_name]
        for _ in cell.domains:
            print(' ', end='')

        print('|', end='')

    print()

    previous_domains = 0
    print('|', end='')
    for cell_name in cells:
        cell = cell_types[cell_name]
        for i in range(len(cell.domains)):
            coverings = register.get_coverings_at_domain_index(previous_domains + i)
            if len(coverings) == 0:
                print(' ', end='')
            elif len(coverings) == 1:
                strand = strand_types[coverings[0]['strand_name']]
                index = previous_domains + i - coverings[0]['start_index']
                if cell.domains[i] == strand.domains[index]:
                    print(' ', end='')
                else:
                    print('↗', end='')
            else:
                print('↗', end='')

        print('|', end='')
        previous_domains += len(cell.domains)

    if len(coverings) > 1:
        last_covering = coverings[-1]
        strand = strand_types[last_covering['strand_name']]
        for _ in range(previous_domains, last_covering['start_index'] + len(strand.domains)):
            print('↗', end='')

    print()

    previous_domains = 0
    print('|', end='')
    for cell_name in cells:
        cell = cell_types[cell_name]
        for i in range(len(cell.domains)):
            coverings = register.get_coverings_at_domain_index(previous_domains + i)
            if len(coverings) == 0:
                print('□', end='')
            elif len(coverings) == 1:
                strand = strand_types[coverings[0]['strand_name']]
                index = previous_domains + i - coverings[0]['start_index']
                if cell.domains[i] == strand.domains[index]:
                    if index < len(strand.domains) - 1:
                        print('=', end='')
                    else:
                        print('>', end='')
                else:
                    print('⭜', end='')
            else:
                print('⭜', end='')

        print('|', end='')
        previous_domains += len(cell.domains)

    print()


def save_data():
    filename = input("Enter filename: ")
    with open(filename, 'w') as file:
        json.dump({
            'cell_types': cell_types,
            'strand_types': strand_types,
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
