import sys

from simd_dna import *
import copy
import json

program_loop = True
local_simulation = Simulation()


def add_cell_type():
    name = input('Enter cell type name: ')
    domains = input('Enter domain names, separated by commas: ').split(',')
    local_simulation.add_cell_type(name, domains)


def add_cells_to_register():
    register_name = input('Enter register name: ')
    cell_name = input('Enter the cell name: ')
    if cell_name not in local_simulation.cell_types.keys():
        print('No such cell exists')
        return

    copies = int(input('Enter the number of copies of each cell: '))
    coverings = input('Enter initial coverings for each cell copy: ')

    if coverings == '':
        coverings = []
    else:
        coverings = coverings.split(';')

    coverings = list(map(extract_covering_offset_pair, coverings))

    local_simulation.add_cells_to_register(register_name, cell_name, coverings, copies)


def extract_covering_offset_pair(value):
    pair = value.split(',')
    pair[1] = int(pair[1])
    return pair


def add_strand_type():
    name = input('What is the strand name? ')
    domains = input('Enter domains separated by commas: ').split(',')
    is_complementary = input('Is the strand complementary to the top strand? Y or N? ').lower() == 'y'
    color = input('What is the color hex code? ')
    if color == '':
        local_simulation.add_strand_type(name, domains, is_complementary)
    else:
        local_simulation.add_strand_type(name, domains, is_complementary, color)


def add_instruction():
    instruction = []
    num_instructions = int(input('Enter the number of instruction strands: '))
    for _ in range(num_instructions):
        strand_name = input('Enter the strand name: ')
        instruction.append(strand_name)

    local_simulation.add_instruction(instruction)


def add_cell_strand_label():
    data = input('Enter the cell name, followed by the coordinate-strand name pairs, followed by the string label, \
all separated by commas: ').split(',')
    coordinate_strand_pairs = []
    for i in range(1, len(data) - 1, 2):
        index = int(data[i])
        coordinate_strand_pairs.append([index, data[i + 1]])
    local_simulation.add_cell_strand_label(data[0], coordinate_strand_pairs, data[-1])


def run_simulation():
    for register_key in local_simulation.registers.keys():
        print(register_key)
        register = local_simulation.registers[register_key]
        if local_simulation.keep_results:
            original_register = None
        else:
            original_register = copy.deepcopy(register)
            local_simulation.keep_results = True  # temporarily save results during instruction cycle
        register.svg_initialize(register_key, len(local_simulation.instructions))

        for inst_num in range(len(local_simulation.instructions)):
            label = ("" if Register.compress_svg_drawings else "Instruction ") + str(inst_num + 1)
            register, before_register, new_strands, unattached_matches = local_simulation.run_instruction(register_key,
                                                                                                          inst_num)

            print("Instruction", inst_num + 1)

            if (len(new_strands) == 0 and (
                    unattached_matches is None or len(unattached_matches) == 0)) and inst_num > 0:
                print('No changes\n')
            else:
                before_register.print(new_strands, unattached_matches)
                print()

            if len(new_strands) > 0:
                before_register._dwg = register._dwg
                before_register.svg_draw_contents(label, len(new_strands) == 0)
                register.svg_draw_strands(new_strands, 3)
                register.svg_draw_strands(unattached_matches, 3 if Register.compress_svg_drawings else 6, True)
                register.svg_increment_vertical_offset()

            if local_simulation.step_by_step_simulation:
                input('Press Enter to continue')

        print("Final result")
        register.print()
        print()
        label = "F" if Register.compress_svg_drawings else "Final result"
        register.svg_draw_contents(label=label)
        register.save_svg()
        if original_register is not None:
            local_simulation.registers[register_key] = original_register
            local_simulation.keep_results = False

        if local_simulation.step_by_step_simulation:
            input('Press Enter to continue')


def save_data():
    filename = input("Enter filename: ")
    register_copies = copy.deepcopy(local_simulation.registers)
    for register in register_copies.values():
        del register.cell_types
        del register.strand_types
        del register._dwg
        del register._svg_current_size_parameters
        del register._svg_vertical_offset
        del register._total_domains

    with open(filename, 'w') as file:
        json.dump({
            'cell_types': local_simulation.cell_types,
            'strand_types': local_simulation.strand_types,
            'registers': register_copies,
            'instructions': local_simulation.instructions
        }, file, indent=4, cls=ObjectEncoder)
        file.flush()


def toggle_step_by_step_simulation():
    local_simulation.step_by_step_simulation = not local_simulation.step_by_step_simulation


def toggle_keep_results():
    local_simulation.keep_results = not local_simulation.keep_results


def toggle_show_unused_instruction_strands():
    local_simulation.show_unused_instruction_strands = not local_simulation.show_unused_instruction_strands


def toggle_compress_svg_drawings():
    Register.compress_svg_drawings = not Register.compress_svg_drawings


def convert_tm_to_simd_wrapper():
    convert_tm_to_simd(local_simulation)


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
            local_simulation.cell_types = decode_json_dict(data['cell_types'], Cell)
            local_simulation.strand_types = decode_json_dict(data['strand_types'], Strand)
            for key in data['registers']:
                data['registers'][key]['cell_types'] = local_simulation.cell_types
                data['registers'][key]['strand_types'] = local_simulation.strand_types
            local_simulation.registers = decode_json_dict(data['registers'], Register)
            local_simulation.instructions = data['instructions']

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
                   '12': convert_tm_to_simd_wrapper,
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
8 - Turn step-by-step simulation ''' + ('off\n' if local_simulation.step_by_step_simulation else 'on\n') +
                       '''9 - ''' + (
                           'Don\'t keep results after simulation\n' if local_simulation.keep_results else 'Keep results after simulation\n') +
                       '''10 - ''' + ('Don\'t show unused instruction strands\n'
                                      if local_simulation.show_unused_instruction_strands
                                      else 'Show unused instruction strands\n') +
                       '''11 - ''' + ('Don\'t compress SVG drawings\n' if Register.compress_svg_drawings
                                      else 'Compress SVG drawings\n') +
                       '''12 - Convert turingmachine.io Turing machine to SIMD register
13 - Exit

''')

        if choice in choice_dict:
            choice_dict[choice]()


if __name__ == '__main__':
    simd_simulator(sys.argv)
