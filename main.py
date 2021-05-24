from simd_dna import *
import copy
import sys
import json
import re
import ruamel.yaml as yaml

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


def convert_tm_to_simd():
    filename = input('Enter the YAML file name: ')
    with open(filename, 'r') as file:
        data = yaml.safe_load(file)
        if not 'input' in data.keys():
            data['input'] = ''

        blank = str(data['blank'])

        state_table_contains_unsupported_character = False
        transition_data = {}
        table = data['table']
        for state in table.keys():
            if state_table_contains_unsupported_character:
                break

            transition = table[state]
            if transition is not None:
                for key in transition.keys():
                    if 'write' in transition[key].keys() \
                            and str(transition[key]['write']) not in ['0', '1', blank]:
                        state_table_contains_unsupported_character = True
                        break

                    if type(key) is tuple:
                        if contains_outside_list(key, ['0', '1', blank]):
                            state_table_contains_unsupported_character = True
                            break
                        else:
                            for item in key:
                                add_simd_transition(state, str(item), copy.deepcopy(transition[key]), transition_data)
                    else:
                        if str(key) not in ['0', '1', blank]:
                            state_table_contains_unsupported_character = True
                            break
                        else:
                            add_simd_transition(state, str(key), transition[key], transition_data)

        if re.sub(rf'[{blank}01]', '', data['input']) != '' or state_table_contains_unsupported_character:
            print('Sorry, the TM conversion algorithm only supports the blank, 0, and 1 alphabet symbols.')
            return

        register_name = input('What name would you like to give the register? ')
        generate_tm_to_simd_data_from_transitions(transition_data, data, register_name)


def contains_outside_list(items, restricted_chars):
    for item in items:
        if str(item) not in restricted_chars:
            return True

    return False


def add_simd_transition(state, read, data, transition_data):
    key = (state, read)
    if 'write' not in data.keys():
        data['write'] = read
    else:
        data['write'] = str(data['write'])

    transition_data[key] = data


def generate_tm_to_simd_data_from_transitions(transition_data, tm_data, register_name):
    global instructions
    Cell.cell_types = {}
    Strand.strand_types = {}
    Register.registers = {}
    instructions = []

    # Create cell type
    tape_cell = Cell([])
    domain_template = '({},{})-{}'
    for key in transition_data.keys():
        for i in range(1, 3):
            tape_cell.domains.append(domain_template.format(*key, i))

    for i in range(1, 9):
        tape_cell.domains.append(str(i))

    Cell.cell_types['Tape cell'] = tape_cell

    create_basic_strand_types(transition_data, domain_template)
    create_tm_cell_labels(transition_data, tm_data['blank'])
    encode_register_data(register_name, transition_data, tm_data)
    generate_tm_instructions(transition_data, tm_data['blank'])


def create_basic_strand_types(transition_data, domain_template):
    import seaborn as sns
    palette = sns.color_palette(None, len(transition_data) + 3)
    palette = [convert_rgb_to_hex(*x) for x in palette]

    strand_name_template = '({},{})_{}'
    strand_data = Strand([], False)
    for key, color in zip(transition_data.keys(), palette[:-3]):
        full_key = strand_name_template.format(*key, 'full')
        Strand.strand_types[full_key] = copy.deepcopy(strand_data)
        for i in range(1, 3):
            Strand.strand_types[full_key].domains.append(domain_template.format(*key, i))
            Strand.strand_types[full_key].color = color

    # Strand patterns
    # Blank - 123,456,78
    # Zero - 123,45,678
    # One - 12,345678

    blank_color = palette[-3]
    zero_color = palette[-2]
    one_color = palette[-1]
    Strand.strand_types['symbol_covered'] = copy.deepcopy(strand_data)
    for i in range(1, 9):
        Strand.strand_types['symbol_covered'].domains.append(str(i))

    strand_data.color = one_color
    Strand.strand_types['symbol_12'] = copy.deepcopy(strand_data)
    for i in range(1, 3):
        Strand.strand_types['symbol_12'].domains.append(str(i))

    Strand.strand_types['symbol_345678'] = copy.deepcopy(strand_data)
    for i in range(3, 9):
        Strand.strand_types['symbol_345678'].domains.append(str(i))

    strand_data.color = zero_color
    Strand.strand_types['symbol_123'] = copy.deepcopy(strand_data)
    for i in range(1, 4):
        Strand.strand_types['symbol_123'].domains.append(str(i))

    Strand.strand_types['symbol_45'] = copy.deepcopy(strand_data)
    for i in range(4, 6):
        Strand.strand_types['symbol_45'].domains.append(str(i))

    Strand.strand_types['symbol_678'] = copy.deepcopy(strand_data)
    for i in range(6, 9):
        Strand.strand_types['symbol_678'].domains.append(str(i))

    strand_data.color = blank_color
    Strand.strand_types['symbol_456'] = copy.deepcopy(strand_data)
    for i in range(4, 7):
        Strand.strand_types['symbol_456'].domains.append(str(i))

    Strand.strand_types['symbol_78'] = copy.deepcopy(strand_data)
    for i in range(7, 9):
        Strand.strand_types['symbol_78'].domains.append(str(i))


def create_tm_cell_labels(transition_data, blank_symbol):
    label_template = {'strands': [], 'label': ''}
    current_index = 0
    for configuration in transition_data.keys():
        label_template['strands'].append([current_index, '({},{})_full'.format(*configuration)])
        current_index += 2

    blank_label = copy.deepcopy(label_template)
    blank_label['label'] = blank_symbol if blank_symbol != ' ' else '⎵'
    blank_label['strands'].append([current_index, 'symbol_123'])
    blank_label['strands'].append([current_index + 3, 'symbol_456'])
    blank_label['strands'].append([current_index + 6, 'symbol_78'])
    Cell.cell_types['Tape cell'].strand_labels.append(blank_label)

    zero_label = copy.deepcopy(label_template)
    zero_label['label'] = '0'
    zero_label['strands'].append([current_index, 'symbol_123'])
    zero_label['strands'].append([current_index + 3, 'symbol_45'])
    zero_label['strands'].append([current_index + 5, 'symbol_678'])
    Cell.cell_types['Tape cell'].strand_labels.append(zero_label)

    one_label = copy.deepcopy(label_template)
    one_label['label'] = '1'
    one_label['strands'].append([current_index, 'symbol_12'])
    one_label['strands'].append([current_index + 2, 'symbol_345678'])
    Cell.cell_types['Tape cell'].strand_labels.append(one_label)

    label_template['strands'].append([current_index, 'symbol_covered'])
    for i in range(len(transition_data)):
        data = list(transition_data.keys())[i]
        if data[1] == ' ':
            data = (data[0], '⎵')
        label_string = '({},{})'.format(*data)

        open_label = copy.deepcopy(label_template)
        open_label['strands'].pop(i)
        open_label['label'] = label_string
        Cell.cell_types['Tape cell'].strand_labels.append(open_label)

        plug_label = copy.deepcopy(label_template)
        plug_label['strands'][i][1] = plug_label['strands'][i][1].replace('full', 'plug')
        plug_label['label'] = label_string
        Cell.cell_types['Tape cell'].strand_labels.append(plug_label)

        final_label = copy.deepcopy(label_template)
        final_label['strands'][i][0] += -1
        final_label['strands'][i][1] = final_label['strands'][i][1].replace('full', 'final')
        final_label['label'] = label_string
        Cell.cell_types['Tape cell'].strand_labels.append(final_label)


def encode_register_data(register_name, transition_data, tm_data):
    Register.registers[register_name] = Register()
    if len(tm_data['input']) > 0:
        Register.registers[register_name].add_cell('Tape cell')
        initial_symbol = tm_data['input'][0]
        initial_configuration = (tm_data['start state'], initial_symbol)
        current_index = 0
        has_valid_initial_transition = False
        for configuration in transition_data.keys():
            if initial_configuration == configuration:
                has_valid_initial_transition = True
            else:
                Register.registers[register_name].coverings.append({
                    'start_index': current_index,
                    'strand_name': '({},{})_full'.format(*configuration)
                })

            current_index += 2
        if has_valid_initial_transition:
            Register.registers[register_name].coverings.append({
                'start_index': current_index,
                'strand_name': 'symbol_covered'
            })
        elif initial_symbol == tm_data['blank']:
            insert_blank_symbol(Register.registers[register_name].coverings, current_index)
        elif initial_symbol == '0':
            insert_zero_symbol(Register.registers[register_name].coverings, current_index)
        else:
            insert_one_symbol(Register.registers[register_name].coverings, current_index)
        current_index += 8
        for i in range(1, len(tm_data['input'])):
            symbol = tm_data['input'][i]
            Register.registers[register_name].add_cell('Tape cell')
            for configuration in transition_data.keys():
                Register.registers[register_name].coverings.append({
                    'start_index': current_index,
                    'strand_name': '({},{})_full'.format(*configuration)
                })
                current_index += 2
            if symbol == tm_data['blank']:
                insert_blank_symbol(Register.registers[register_name].coverings, current_index)
            elif symbol == '0':
                insert_zero_symbol(Register.registers[register_name].coverings, current_index)
            else:
                insert_one_symbol(Register.registers[register_name].coverings, current_index)
            current_index += 8


def insert_blank_symbol(coverings, current_index):
    coverings.append({
        'start_index': current_index,
        'strand_name': 'symbol_123'
    })
    current_index += 3
    coverings.append({
        'start_index': current_index,
        'strand_name': 'symbol_456'
    })
    current_index += 3
    coverings.append({
        'start_index': current_index,
        'strand_name': 'symbol_78'
    })


def insert_zero_symbol(coverings, current_index):
    coverings.append({
        'start_index': current_index,
        'strand_name': 'symbol_123'
    })
    current_index += 3
    coverings.append({
        'start_index': current_index,
        'strand_name': 'symbol_45'
    })
    current_index += 2
    coverings.append({
        'start_index': current_index,
        'strand_name': 'symbol_678'
    })


def insert_one_symbol(coverings, current_index):
    coverings.append({
        'start_index': current_index,
        'strand_name': 'symbol_12'
    })
    current_index += 2
    coverings.append({
        'start_index': current_index,
        'strand_name': 'symbol_345678'
    })


def generate_tm_instructions(transition_data, blank_symbol):
    global instructions
    configurations = list(transition_data.keys())
    for configuration in configurations:
        domains = []
        color = Strand.strand_types['({},{})_full'.format(*configuration)].color
        for i in range(1, 3):
            domains.append('({},{})-{}'.format(*configuration, i))
        domains.append('9')
        plug = Strand(domains, False, color)
        Strand.strand_types['({},{})_plug'.format(*configuration)] = plug
        unplug = copy.deepcopy(plug)
        unplug.is_complementary = True
        Strand.strand_types['({},{})_unplug'.format(*configuration)] = unplug

        domains = copy.deepcopy(domains)
        domains.pop()
        domains.insert(0, '10')
        final = Strand(domains, False, color)
        Strand.strand_types['({},{})_final'.format(*configuration)] = final
        restart = copy.deepcopy(final)
        restart.is_complementary = True
        Strand.strand_types['({},{})_restart'.format(*configuration)] = restart

    instructions.append(['({},{})_plug'.format(*configuration) for configuration in configurations])

    num_configurations = len(transition_data)
    for i in range(num_configurations):
        configuration = configurations[i]
        instructions.append(['({},{})_unplug'.format(*configuration)])

        transition = transition_data[configuration]
        if 'L' in transition.keys():
            create_left_instruction(i, num_configurations, transition, transition_data, blank_symbol)
        else:
            create_right_instruction(i, num_configurations, transition, transition_data, blank_symbol)

    # Restart instruction
    instructions.append(['({},{})_restart'.format(*c) for c in configurations])


def create_left_instruction(index, num_configurations, transition, transition_data, blank_symbol):
    global instructions
    next_state = transition['L']
    write_symbol = transition['write']
    configurations = list(transition_data.keys())

    # Displace all strands until left symbol area
    if index > 0:
        left_instruction = ['left_start']
        if 'left_start' not in Strand.strand_types.keys():
            first_domain = '({},{})-1'.format(*configurations[0])
            domains = ['7', '8', first_domain, '11']
            Strand.strand_types['left_start'] = Strand(domains, False)
            Strand.strand_types['left_start_strip'] = Strand(domains, True)
    else:
        left_instruction = ['left_start_extended']
        if 'left_start_extended' not in Strand.strand_types.keys():
            domains = ['7', '8']
            for i in range(1, 3):
                domains.append('({},{})-{}'.format(*configurations[0], i))
            domains.append('11')
            Strand.strand_types['left_start_extended'] = Strand(domains, False)
            Strand.strand_types['left_start_extended_strip'] = Strand(domains, True)
    for i in range(0, index - 1):
        current_configuration = configurations[i]
        next_configuration = configurations[i + 1]
        left_strand_name = '({},{})_left'.format(*current_configuration)
        if left_strand_name not in Strand.strand_types.keys():
            left_strand_domains = ['({},{})-2'.format(*current_configuration),
                                   '({},{})-1'.format(*next_configuration), '11']
            color = Strand.strand_types['({},{})_full'.format(*current_configuration)].color
            Strand.strand_types[left_strand_name] = Strand(left_strand_domains, False, color)
        left_instruction.append(left_strand_name)
    if index > 0:
        current_configuration = configurations[index - 1]
        next_configuration = configurations[index]
        left_extended_name = '({},{})_left_extended'.format(*current_configuration)
        if left_extended_name not in Strand.strand_types.keys():
            left_extended_domains = ['({},{})-2'.format(*current_configuration)]
            for i in range(1, 3):
                left_extended_domains.append('({},{})-{}'.format(*next_configuration, i))
            left_extended_domains.append('11')
            color = Strand.strand_types['({},{})_full'.format(*current_configuration)].color
            Strand.strand_types[left_extended_name] = Strand(left_extended_domains, False, color)
            left_extended_strip_name = '({},{})_left_extended_strip'.format(*current_configuration)
            Strand.strand_types[left_extended_strip_name] = Strand(left_extended_domains, True, color)
        left_instruction.append(left_extended_name)
    instructions.append(left_instruction)

    if 'one_left_start' not in Strand.strand_types.keys():
        first_domain = '({},{})-1'.format(*configurations[0])
        color = Strand.strand_types['symbol_345678'].color
        domains = ['3', '4', '5', '6', '7', '8', first_domain, '13']
        Strand.strand_types['one_left_start'] = Strand(domains, False, color)
        Strand.strand_types['one_left_strip'] = Strand(domains, True, color)

    # Add instructions for when left cell is zero
    next_index = configurations.index((next_state, '0')) \
        if (next_state, '0') in configurations else -1
    if next_index == -1:
        if 'temp_zero' not in Strand.strand_types.keys():
            color = Strand.strand_types['symbol_45'].color
            domains = ['4', '5', '6', '13']
            Strand.strand_types['temp_zero'] = Strand(domains, False, color)
            Strand.strand_types['temp_zero_strip'] = Strand(domains, True, color)
        instructions.append(['temp_zero'])
        instructions.append(['one_left_start'])
        instructions.append(['temp_zero_strip'])
        instructions.append(['symbol_45', 'symbol_678'])
    else:
        zero_instruction = []
        generate_left_cascade_instruction_strands(next_index, num_configurations, configurations, zero_instruction)
        if 'zero_symbol_displace' not in Strand.strand_types.keys():
            domains = ['2', '3', '4', '5', '6', '12']
            color = Strand.strand_types['symbol_45'].color
            Strand.strand_types['zero_symbol_displace'] = Strand(domains, False, color)
        zero_instruction.append('zero_symbol_displace')
        instructions.append(zero_instruction)
        instructions.append(['one_left_start'])  # in case left cell is 1
        generate_left_final_instruction_strands(next_index, num_configurations, configurations)

    # Add instructions for when left cell is one
    next_index = configurations.index((next_state, '1')) \
        if (next_state, '1') in configurations else -1
    instructions.append(['one_left_strip'])
    if next_index == -1:
        instructions.append(['symbol_345678'])
    else:
        if 'temp_one' not in Strand.strand_types.keys():
            color = Strand.strand_types['symbol_345678'].color
            first_domain = '({},{})-1'.format(*configurations[0])
            domains = ['4', '5', '6', '7', '8', first_domain, '13']
            Strand.strand_types['temp_one'] = Strand(domains, False, color)

        instructions.append(['temp_one'])

        one_instruction = []
        generate_left_cascade_instruction_strands(next_index, num_configurations, configurations, one_instruction)
        if 'one_symbol_displace' not in Strand.strand_types.keys():
            domains = ['2', '3', '12']
            color = Strand.strand_types['symbol_12'].color
            Strand.strand_types['one_symbol_displace'] = Strand(domains, False, color)
        one_instruction.append('one_symbol_displace')
        instructions.append(one_instruction)
        generate_left_final_instruction_strands(next_index, num_configurations, configurations)

    # Add instructions for when left cell is blank
    next_index = configurations.index((next_state, blank_symbol)) \
        if (next_state, blank_symbol) in configurations else -1
    if index > 0:
        instructions.append(['left_start_strip'])
    else:
        instructions.append(['left_start_extended_strip'])

    if next_index == -1:
        instructions.append(['symbol_78'])
    else:
        instructions.append(['symbol_right'])
        instructions.append(['symbol_right_strip'])
        blank_instruction = []
        generate_left_cascade_instruction_strands(next_index, num_configurations, configurations, blank_instruction)
        instructions.append(blank_instruction)
        generate_left_final_instruction_strands(next_index, num_configurations, configurations)

    # Replace current cell coverings
    if index > 0:
        instructions.append(['({},{})_left_extended_strip'.format(*configurations[index - 1])])
    else:
        instructions.append(['left_start_extended_strip'])

    replacement_instructions = []
    for i in range(0, index):
        replacement_instructions.append('({},{})_full'.format(*configurations[i]))
    generate_right_cascade_instruction_strands(index + 1, num_configurations, configurations, replacement_instructions)
    instructions.append(replacement_instructions)

    right_strip_instructions = []
    for i in range(index + 1, num_configurations):
        current_configuration = configurations[i]
        right_strip_instructions.append('({},{})_right_strip'.format(*current_configuration))
    right_strip_instructions.append('symbol_right_strip')
    instructions.append(right_strip_instructions)

    last_one_instructions = []
    for i in range(index, num_configurations):
        current_configuration = configurations[i]
        last_one_instructions.append('({},{})_full'.format(*current_configuration))
    if write_symbol == blank_symbol:
        last_one_instructions.extend(['symbol_123', 'symbol_456', 'symbol_78'])
    elif write_symbol == '0':
        last_one_instructions.extend(['symbol_123', 'symbol_45', 'symbol_678'])
    else:
        last_one_instructions.extend(['symbol_12', 'symbol_345678'])

    instructions.append(last_one_instructions)


def generate_left_cascade_instruction_strands(next_index, num_configurations, configurations, instruction_list):
    for i in range(next_index, num_configurations):
        current_configuration = configurations[i]
        next_configuration = configurations[i + 1] if i < num_configurations - 1 else -1
        left_strand_name = '({},{})_left'.format(*current_configuration)
        if left_strand_name not in Strand.strand_types.keys():
            left_strand_domains = ['({},{})-2'.format(*current_configuration)]
            if next_configuration == -1:
                left_strand_domains.append('1')
            else:
                left_strand_domains.append('({},{})-{}'.format(*next_configuration, 1))
            left_strand_domains.append('11')
            color = Strand.strand_types['({},{})_full'.format(*current_configuration)].color
            Strand.strand_types[left_strand_name] = Strand(left_strand_domains, False, color)
        instruction_list.append(left_strand_name)


def generate_left_final_instruction_strands(next_index, num_configurations, configurations):
    global instructions
    new_cell_instructions = ['({},{})_final'.format(*configurations[next_index])]
    for i in range(next_index + 1, num_configurations):
        new_cell_instructions.append('({},{})_full'.format(*configurations[i]))
    new_cell_instructions.append('symbol_covered')
    instructions.append(new_cell_instructions)


def create_right_instruction(index, num_configurations, transition, transition_data, blank_symbol):
    global instructions
    next_state = transition['R']
    write_symbol = transition['write']
    configurations = list(transition_data.keys())

    # Displace all strands until right end of cell
    right_instruction = []

    generate_right_cascade_instruction_strands(index + 1, num_configurations, configurations, right_instruction)
    instructions.append(right_instruction)

    # Strip out strands from previous instruction
    right_strip_instructions = []
    for i in range(index + 1, num_configurations):
        current_configuration = configurations[i]
        right_strip_instructions.append('({},{})_right_strip'.format(*current_configuration))
    right_strip_instructions.append('symbol_right_strip')

    instructions.append(right_strip_instructions)

    # Replace current cell coverings
    replacement_instructions = []
    for current_configuration in configurations:
        replacement_instructions.append('({},{})_full'.format(*current_configuration))

    if write_symbol == blank_symbol:
        replacement_instructions.append('symbol_123')
        if 'symbol_45(13)' not in Strand.strand_types.keys():
            domains = ['4', '5', '13']
            Strand.strand_types['symbol_45(13)'] = Strand(domains, False, Strand.strand_types['symbol_456'].color)
            Strand.strand_types['symbol_45(13)_strip'] = Strand(domains, True, Strand.strand_types['symbol_456'].color)
            domains = ['6', '7', '13']
            Strand.strand_types['symbol_67(13)'] = Strand(domains, False, Strand.strand_types['symbol_78'].color)
            Strand.strand_types['symbol_67(13)_strip'] = Strand(domains, True, Strand.strand_types['symbol_78'].color)
        replacement_instructions.extend(['symbol_45(13)', 'symbol_67(13)'])
    elif write_symbol == '0':
        replacement_instructions.append('symbol_123')
        replacement_instructions.append('symbol_45')
        if 'symbol_67' not in Strand.strand_types.keys():
            Strand.strand_types['symbol_67'] = Strand(['6', '7'], False, Strand.strand_types['symbol_678'].color)
        replacement_instructions.append('symbol_67')
    else:
        replacement_instructions.append('symbol_12')
        if 'symbol_34567' not in Strand.strand_types.keys():
            Strand.strand_types['symbol_34567'] = Strand(['3', '4', '5', '6', '7'], False,
                                                         Strand.strand_types['symbol_345678'].color)
        replacement_instructions.append('symbol_34567')

    instructions.append(replacement_instructions)

    # Displace right cell contents
    next_right_instructions = []
    if '({},{})_next_right'.format(*configurations[0]) not in Strand.strand_types.keys():
        for i in range(len(configurations)):
            domains = ['8' if i == 0 else '({},{})-2'.format(*configurations[i - 1])]
            current_configuration = configurations[i]
            domains.append('({},{})-1'.format(*current_configuration))
            domains.append('11')
            strand_name = '({},{})_next_right'.format(*current_configuration)
            color = Strand.strand_types['({},{})_full'.format(*current_configuration)].color
            Strand.strand_types[strand_name] = Strand(domains, False, color)
            strand_strip_name = '({},{})_next_right_strip'.format(*current_configuration)
            Strand.strand_types[strand_strip_name] = Strand(domains, True, color)

        last_configuration = configurations[num_configurations - 1]
        domains = ['({},{})-2'.format(*last_configuration), '1', '2', '11']
        Strand.strand_types['next_right_1'] = Strand(domains, False)
        Strand.strand_types['next_right_1_strip'] = Strand(domains, True)
        Strand.strand_types['next_right_2'] = Strand(['3', '4', '5', '11'], False)
        Strand.strand_types['next_right_2_strip'] = Strand(Strand.strand_types['next_right_2'].domains, True)

    for current_configuration in configurations:
        next_right_instructions.append('({},{})_next_right'.format(*current_configuration))
    next_right_instructions.append('next_right_1')
    next_right_instructions.append('next_right_2')

    instructions.append(next_right_instructions)

    # Add instructions for when right cell is blank
    next_index = configurations.index((next_state, blank_symbol)) \
        if (next_state, blank_symbol) in configurations else -1
    first_configuration = configurations[0]
    if 'right_blank' not in Strand.strand_types.keys():
        domains = ['4', '5', '6', '12']
        Strand.strand_types['right_blank'] = Strand(domains, False)
        Strand.strand_types['right_blank_strip'] = Strand(domains, True)
        domains = ['8']
        for i in range(1, 3):
            domains.append('({},{})-{}'.format(*first_configuration, i))
        domains.append('12')
        color = Strand.strand_types['({},{})_full'.format(*first_configuration)].color
        Strand.strand_types['8_right_plug'] = Strand(domains, False, color)
        Strand.strand_types['8_right_plug_strip'] = Strand(domains, True, color)
    blank_instructions = ['8_right_plug']
    for i in range(1, len(configurations)):
        current_configuration = configurations[i]
        strand_type = '({},{})_full' if i != next_index else '({},{})_final'
        blank_instructions.append(strand_type.format(*current_configuration))
    blank_instructions.append('symbol_123')
    blank_instructions.append('right_blank')
    instructions.append(blank_instructions)
    instructions.append(['right_blank_strip'])
    instructions.append(['symbol_456' if next_index == -1 else 'symbol_covered'])
    instructions.append(['8_right_plug_strip'])
    last_blank_instruction = []
    if write_symbol == blank_symbol:
        last_blank_instruction.extend(['symbol_456', 'symbol_78'])
    elif write_symbol == '0':
        last_blank_instruction.append('symbol_678')
    else:
        last_blank_instruction.append('symbol_345678')
    strand_type = '({},{})_final' if next_index == 0 else '({},{})_full'
    last_blank_instruction.append(strand_type.format(*first_configuration))
    instructions.append(last_blank_instruction)

    # Add instructions for when right cell is zero
    next_index = configurations.index((next_state, '0')) \
        if (next_state, '0') in configurations else -1
    instructions.append(['next_right_2_strip'])
    if next_index == -1:
        instructions.append(['symbol_123', 'symbol_45'])
    else:
        instructions.append(['symbol_covered'])
    zero_instructions = ['8_right_plug']
    for i in range(1, len(configurations)):
        current_configuration = configurations[i]
        strand_type = '({},{})_full' if i != next_index else '({},{})_final'
        zero_instructions.append(strand_type.format(*current_configuration))
    instructions.append(zero_instructions)
    instructions.append(['8_right_plug_strip'])
    last_zero_instruction = []
    if write_symbol == blank_symbol:
        last_zero_instruction.extend(['symbol_456', 'symbol_78'])
    elif write_symbol == '0':
        last_zero_instruction.append('symbol_678')
    else:
        last_zero_instruction.append('symbol_345678')
    strand_type = '({},{})_final' if next_index == 0 else '({},{})_full'
    last_zero_instruction.append(strand_type.format(*first_configuration))
    instructions.append(last_zero_instruction)

    # Add instructions for when right cell is one
    next_index = configurations.index((next_state, '1')) \
        if (next_state, '1') in configurations else -1
    one_instructions = []
    for current_configuration in configurations:
        one_instructions.append('({},{})_next_right_strip'.format(*current_configuration))
    one_instructions.append('next_right_1_strip')
    instructions.append(one_instructions)
    next_one_instruction = []
    if write_symbol == blank_symbol:
        next_one_instruction.extend(['symbol_456', 'symbol_78'])
    elif write_symbol == '0':
        next_one_instruction.append('symbol_678')
    else:
        next_one_instruction.append('symbol_345678')
    for i in range(len(configurations)):
        current_configuration = configurations[i]
        strand_type = '({},{})_full' if i != next_index else '({},{})_final'
        next_one_instruction.append(strand_type.format(*current_configuration))
    instructions.append(next_one_instruction)

    if next_index == -1:
        instructions.append(['symbol_12', 'symbol_345678'])
    else:
        instructions.append(['symbol_covered'])


def generate_right_cascade_instruction_strands(index, num_configurations, configurations, instruction_list):
    for i in range(index, num_configurations):
        previous_configuration = configurations[i - 1] if i > 0 else -1
        current_configuration = configurations[i]
        right_strand_name = '({},{})_right'.format(*current_configuration)
        if right_strand_name not in Strand.strand_types.keys():
            if previous_configuration == -1:
                right_strand_domains = ['8']
            else:
                right_strand_domains = ['({},{})-2'.format(*previous_configuration)]
            right_strand_domains.append('({},{})-1'.format(*current_configuration))
            right_strand_domains.append('11')
            color = Strand.strand_types['({},{})_full'.format(*current_configuration)].color
            Strand.strand_types[right_strand_name] = Strand(right_strand_domains, False, color)
            right_strip_name = '({},{})_right_strip'.format(*current_configuration)
            Strand.strand_types[right_strip_name] = Strand(right_strand_domains, True, color)
        instruction_list.append(right_strand_name)

    if 'symbol_right' not in Strand.strand_types.keys():
        last_configuration = configurations[num_configurations - 1]
        right_strand_domains = ['({},{})-2'.format(*last_configuration)]
        for i in range(1, 8):
            right_strand_domains.append(str(i))
        right_strand_domains.append('11')
        color = Strand.strand_types['symbol_covered'].color
        Strand.strand_types['symbol_right'] = Strand(right_strand_domains, False, color)
        Strand.strand_types['symbol_right_strip'] = Strand(right_strand_domains, True, color)

    instruction_list.append('symbol_right')


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
                   '12': convert_tm_to_simd,
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
