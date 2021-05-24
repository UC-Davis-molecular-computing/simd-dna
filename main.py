from simd_dna import *
import copy
import sys
import json

program_loop = True
step_by_step_simulation = False
keep_results = False
show_unused_instruction_strands = False
instructions = []


def add_cell_type():
    name = input('Enter cell type name: ')
    domains = input('Enter domain names, separated by commas: ').split(',')
    if name not in Cell.cell_types.keys():
        Cell.cell_types[name] = Cell(domains)


def add_cells_to_register():
    register_name = input('Enter register name: ')
    if register_name not in Register.registers:
        Register.registers[register_name] = Register()

    current_register = Register.registers[register_name]

    cell_name = input('Enter the cell name: ')
    if cell_name not in Cell.cell_types.keys():
        print('No such cell exists')
        return

    copies = int(input('Enter the number of copies of each cell: '))
    coverings = input('Enter initial coverings for each cell copy: ')

    if coverings == '':
        coverings = []
    else:
        coverings = coverings.split(';')

    cell_size = len(Cell.cell_types[cell_name].domains)
    for _ in range(copies):
        current_register.add_cell(cell_name)
        for covering in coverings:
            data = covering.split(',')
            strand_type = data[0]
            offset = int(data[1])
            current_register.attempt_attachment(-cell_size + offset, strand_type)


def add_strand_type():
    name = input('What is the strand name? ')
    domains = input('Enter domains separated by commas: ').split(',')
    is_complementary = input('Is the strand complementary to the top strand? Y or N? ')
    color = input('What is the color hex code? ')
    if color == '':
        color = '#000000'

    Strand.strand_types[name] = Strand(domains, True if is_complementary.lower() == 'y' else False, color)


def add_instruction():
    instruction = []
    num_instructions = int(input('Enter the number of instruction strands: '))
    for _ in range(num_instructions):
        strand_name = input('Enter the strand name: ')
        instruction.append(strand_name)

    instructions.append(instruction)


def add_cell_strand_label():
    data = input('Enter the cell name, followed by the coordinate-strand name pairs, followed by the string label, \
all separated by commas: ').split(',')
    label = {'strands': [], 'label': data[-1]}
    for i in range(1, len(data) - 1, 2):
        strands = label['strands']
        index = int(data[i])
        if len(strands) == 0:
            strands.append([index, data[i + 1]])
        else:
            for j in range(len(strands)):
                if strands[j][0] > index:
                    strands.insert(j, [index, data[i + 1]])
                    break

                if j == len(strands) - 1:
                    strands.append([index, data[i + 1]])
    Cell.cell_types[data[0]].strand_labels.append(label)


def run_simulation():
    register_copies = copy.deepcopy(Register.registers)
    for register_key in register_copies.keys():
        print(register_key)
        register = register_copies[register_key]
        register.svg_initialize(register_key, len(instructions))

        total_domains = 0
        for cell_name in register.cells:
            total_domains += len(Cell.cell_types[cell_name].domains)

        for inst_num in range(len(instructions)):
            if show_unused_instruction_strands:
                unattached_matches = []
            else:
                unattached_matches = None

            pre_instruction_register = copy.deepcopy(register)
            label = ("" if Register.compress_svg_drawings else "Instruction ") + str(inst_num + 1)
            inst = instructions[inst_num]
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
                        new_strands.extend([strand for strand in new_attachments if strand not in displaced_strands])

                    if unattached_matches is not None:
                        unattached_matches.extend(displaced_strands)

            print("Instruction", inst_num + 1)
            new_strands.sort(key=lambda x: x['start_index'])
            if unattached_matches is not None:
                unattached_matches = [strand for strand in unattached_matches if strand not in new_strands]
                unattached_matches = register.sanitize_unattached_strands(unattached_matches, new_strands)

            if (len(new_strands) == 0 and (
                    unattached_matches is None or len(unattached_matches) == 0)) and inst_num > 0:
                print('No changes\n')
            else:
                pre_instruction_register.print(new_strands, unattached_matches)
                print()

            if len(new_strands) > 0:
                pre_instruction_register._dwg = register._dwg
                pre_instruction_register.svg_draw_contents(label)
                register.svg_draw_strands(new_strands, 3)
                register.svg_draw_strands(unattached_matches, 3 if Register.compress_svg_drawings else 6, True)
                register.svg_increment_vertical_offset()

            if step_by_step_simulation:
                input('Press Enter to continue')

        print("Final result")
        register.print()
        print()
        label = "F" if Register.compress_svg_drawings else "Final result"
        register.svg_draw_contents(label=label)
        register.save_svg()
        if step_by_step_simulation:
            input('Press Enter to continue')

        if keep_results:
            Register.registers[register_key] = register


def save_data():
    filename = input("Enter filename: ")
    register_copies = copy.deepcopy(Register.registers)
    for register in register_copies.values():
        del register._dwg
        del register._svg_current_size_parameters
        del register._svg_vertical_offset
        del register._total_domains

    with open(filename, 'w') as file:
        json.dump({
            'cell_types': Cell.cell_types,
            'strand_types': Strand.strand_types,
            'registers': register_copies,
            'instructions': instructions
        }, file, indent=4, cls=ObjectEncoder)
        file.flush()


def toggle_step_by_step_simulation():
    global step_by_step_simulation
    step_by_step_simulation = not step_by_step_simulation


def toggle_keep_results():
    global keep_results
    keep_results = not keep_results


def toggle_show_unused_instruction_strands():
    global show_unused_instruction_strands
    show_unused_instruction_strands = not show_unused_instruction_strands


def toggle_compress_svg_drawings():
    Register.compress_svg_drawings = not Register.compress_svg_drawings


def convert_tm_to_simd_wrapper():
    convert_tm_to_simd(instructions)


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
            global instructions
            Cell.cell_types = decode_json_dict(data['cell_types'], Cell)
            Strand.strand_types = decode_json_dict(data['strand_types'], Strand)
            Register.registers = decode_json_dict(data['registers'], Register)
            instructions = data['instructions']

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
8 - Turn step-by-step simulation ''' + ('off\n' if step_by_step_simulation else 'on\n') +
                       '''9 - ''' + (
                           'Don\'t keep results after simulation\n' if keep_results else 'Keep results after simulation\n') +
                       '''10 - ''' + ('Don\'t show unused instruction strands\n' if show_unused_instruction_strands
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
