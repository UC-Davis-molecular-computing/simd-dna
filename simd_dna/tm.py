from typing import Dict, List

from simd_dna.simulation import *
from simd_dna.functions import convert_rgb_to_hex
from simd_dna.classes import TopStrand

import copy
import re
import ruamel.yaml as yaml


def convert_tm_to_simd(simulation: Simulation) -> None:
    filename = input('Enter the YAML file name: ')
    with open(filename, 'r') as file:
        data = yaml.safe_load(file)
        if 'input' not in data.keys():
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
        generate_tm_to_simd_data_from_transitions(transition_data, data, register_name, simulation)


def contains_outside_list(items: List, restricted_chars: List[str]) -> bool:
    for item in items:
        if str(item) not in restricted_chars:
            return True

    return False


def add_simd_transition(state, read, data, transition_data) -> None:
    key = (state, read)
    if 'write' not in data.keys():
        data['write'] = read
    else:
        data['write'] = str(data['write'])

    transition_data[key] = data


def generate_tm_to_simd_data_from_transitions(transition_data, tm_data, register_name, simulation):
    simulation.cell_types.clear()
    simulation.strand_types.clear()
    simulation.registers.clear()
    simulation.instructions.clear()

    # Create cell type
    tape_cell = Cell([])
    domain_template = '({},{})-{}'
    for key in transition_data.keys():
        for i in range(1, 3):
            tape_cell.domains.append(domain_template.format(*key, i))

    for i in range(1, 9):
        tape_cell.domains.append(str(i))

    simulation.cell_types['Tape cell'] = tape_cell

    create_basic_strand_types(transition_data, domain_template, simulation)
    create_tm_cell_labels(transition_data, tm_data['blank'], simulation)
    encode_register_data(register_name, transition_data, tm_data, simulation)
    generate_tm_instructions(transition_data, tm_data['blank'], simulation)


def create_basic_strand_types(transition_data, domain_template, simulation):
    import seaborn as sns
    palette = sns.color_palette(None, len(transition_data) + 3)
    palette = [convert_rgb_to_hex(*x) for x in palette]

    strand_name_template = '({},{})_{}'
    strand_data = Strand([], False)
    for key, color in zip(transition_data.keys(), palette[:-3]):
        full_key = strand_name_template.format(*key, 'full')
        simulation.strand_types[full_key] = copy.deepcopy(strand_data)
        for i in range(1, 3):
            simulation.strand_types[full_key].domains.append(domain_template.format(*key, i))
            simulation.strand_types[full_key].color = color

    # Strand patterns
    # Blank - 123,456,78
    # Zero - 123,45,678
    # One - 12,345678

    blank_color = palette[-3]
    zero_color = palette[-2]
    one_color = palette[-1]
    simulation.strand_types['symbol_covered'] = copy.deepcopy(strand_data)
    for i in range(1, 9):
        simulation.strand_types['symbol_covered'].domains.append(str(i))

    strand_data.color = one_color
    simulation.strand_types['symbol_12'] = copy.deepcopy(strand_data)
    for i in range(1, 3):
        simulation.strand_types['symbol_12'].domains.append(str(i))

    simulation.strand_types['symbol_345678'] = copy.deepcopy(strand_data)
    for i in range(3, 9):
        simulation.strand_types['symbol_345678'].domains.append(str(i))

    strand_data.color = zero_color
    simulation.strand_types['symbol_123'] = copy.deepcopy(strand_data)
    for i in range(1, 4):
        simulation.strand_types['symbol_123'].domains.append(str(i))

    simulation.strand_types['symbol_45'] = copy.deepcopy(strand_data)
    for i in range(4, 6):
        simulation.strand_types['symbol_45'].domains.append(str(i))

    simulation.strand_types['symbol_678'] = copy.deepcopy(strand_data)
    for i in range(6, 9):
        simulation.strand_types['symbol_678'].domains.append(str(i))

    strand_data.color = blank_color
    simulation.strand_types['symbol_456'] = copy.deepcopy(strand_data)
    for i in range(4, 7):
        simulation.strand_types['symbol_456'].domains.append(str(i))

    simulation.strand_types['symbol_78'] = copy.deepcopy(strand_data)
    for i in range(7, 9):
        simulation.strand_types['symbol_78'].domains.append(str(i))


def create_tm_cell_labels(transition_data, blank_symbol, simulation):
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
    simulation.cell_types['Tape cell'].strand_labels.append(blank_label)

    zero_label = copy.deepcopy(label_template)
    zero_label['label'] = '0'
    zero_label['strands'].append([current_index, 'symbol_123'])
    zero_label['strands'].append([current_index + 3, 'symbol_45'])
    zero_label['strands'].append([current_index + 5, 'symbol_678'])
    simulation.cell_types['Tape cell'].strand_labels.append(zero_label)

    one_label = copy.deepcopy(label_template)
    one_label['label'] = '1'
    one_label['strands'].append([current_index, 'symbol_12'])
    one_label['strands'].append([current_index + 2, 'symbol_345678'])
    simulation.cell_types['Tape cell'].strand_labels.append(one_label)

    label_template['strands'].append([current_index, 'symbol_covered'])
    for i in range(len(transition_data)):
        data = list(transition_data.keys())[i]
        if data[1] == ' ':
            data = (data[0], '⎵')
        label_string = '({},{})'.format(*data)

        open_label = copy.deepcopy(label_template)
        open_label['strands'].pop(i)
        open_label['label'] = label_string
        simulation.cell_types['Tape cell'].strand_labels.append(open_label)

        plug_label = copy.deepcopy(label_template)
        plug_label['strands'][i][1] = plug_label['strands'][i][1].replace('full', 'plug')
        plug_label['label'] = label_string
        simulation.cell_types['Tape cell'].strand_labels.append(plug_label)

        final_label = copy.deepcopy(label_template)
        final_label['strands'][i][0] += -1
        final_label['strands'][i][1] = final_label['strands'][i][1].replace('full', 'final')
        final_label['label'] = label_string
        simulation.cell_types['Tape cell'].strand_labels.append(final_label)


def encode_register_data(register_name: str, transition_data, tm_data: Dict, simulation: Simulation) -> None:
    simulation.registers[register_name] = Register(simulation.cell_types, simulation.strand_types)
    if len(tm_data['input']) > 0:
        simulation.registers[register_name].add_cell('Tape cell')
        initial_symbol = tm_data['input'][0]
        initial_configuration = (tm_data['start state'], initial_symbol)
        current_index = 0
        has_valid_initial_transition = False
        for configuration in transition_data.keys():
            if initial_configuration == configuration:
                has_valid_initial_transition = True
            else:
                simulation.registers[register_name].top_strands.append(TopStrand(current_index,
                                                                                 '({},{})_full'.format(*configuration)))

            current_index += 2
        if has_valid_initial_transition:
            simulation.registers[register_name].top_strands.append(TopStrand(current_index, 'symbol_covered'))
        elif initial_symbol == tm_data['blank']:
            insert_blank_symbol(simulation.registers[register_name].top_strands, current_index)
        elif initial_symbol == '0':
            insert_zero_symbol(simulation.registers[register_name].top_strands, current_index)
        else:
            insert_one_symbol(simulation.registers[register_name].top_strands, current_index)
        current_index += 8
        for i in range(1, len(tm_data['input'])):
            symbol = tm_data['input'][i]
            simulation.registers[register_name].add_cell('Tape cell')
            for configuration in transition_data.keys():
                simulation.registers[register_name].top_strands.append(TopStrand(current_index,
                                                                                 '({},{})_full'.format(*configuration)))
                current_index += 2
            if symbol == tm_data['blank']:
                insert_blank_symbol(simulation.registers[register_name].top_strands, current_index)
            elif symbol == '0':
                insert_zero_symbol(simulation.registers[register_name].top_strands, current_index)
            else:
                insert_one_symbol(simulation.registers[register_name].top_strands, current_index)
            current_index += 8


def insert_blank_symbol(top_strands, current_index):
    top_strands.append(TopStrand(current_index, 'symbol_123'))
    current_index += 3
    top_strands.append(TopStrand(current_index, 'symbol_456'))
    current_index += 3
    top_strands.append(TopStrand(current_index, 'symbol_78'))


def insert_zero_symbol(top_strands, current_index):
    top_strands.append(TopStrand(current_index, 'symbol_123'))
    current_index += 3
    top_strands.append(TopStrand(current_index, 'symbol_45'))
    current_index += 2
    top_strands.append(TopStrand(current_index, 'symbol_678'))


def insert_one_symbol(top_strands, current_index):
    top_strands.append(TopStrand(current_index, 'symbol_12'))
    current_index += 2
    top_strands.append(TopStrand(current_index, 'symbol_345678'))


def generate_tm_instructions(transition_data, blank_symbol, simulation):
    configurations = list(transition_data.keys())
    for configuration in configurations:
        domains = []
        color = simulation.strand_types['({},{})_full'.format(*configuration)].color
        for i in range(1, 3):
            domains.append('({},{})-{}'.format(*configuration, i))
        domains.append('9')
        plug = Strand(domains, False, color)
        simulation.strand_types['({},{})_plug'.format(*configuration)] = plug
        unplug = copy.deepcopy(plug)
        unplug.is_complementary = True
        simulation.strand_types['({},{})_unplug'.format(*configuration)] = unplug

        domains = copy.deepcopy(domains)
        domains.pop()
        domains.insert(0, '10')
        final = Strand(domains, False, color)
        simulation.strand_types['({},{})_final'.format(*configuration)] = final
        restart = copy.deepcopy(final)
        restart.is_complementary = True
        simulation.strand_types['({},{})_restart'.format(*configuration)] = restart

    simulation.instructions.append(['({},{})_plug'.format(*configuration) for configuration in configurations])

    num_configurations = len(transition_data)
    for i in range(num_configurations):
        configuration = configurations[i]
        simulation.instructions.append(['({},{})_unplug'.format(*configuration)])

        transition = transition_data[configuration]
        if 'L' in transition.keys():
            create_left_instruction(i, num_configurations, transition, transition_data, blank_symbol, simulation)
        else:
            create_right_instruction(i, num_configurations, transition, transition_data, blank_symbol, simulation)

    # Restart instruction
    simulation.instructions.append(['({},{})_restart'.format(*c) for c in configurations])


def create_left_instruction(index, num_configurations, transition, transition_data, blank_symbol, simulation):
    next_state = transition['L']
    write_symbol = transition['write']
    configurations = list(transition_data.keys())

    # Displace all strands until left symbol area
    if index > 0:
        left_instruction = ['left_start']
        if 'left_start' not in simulation.strand_types.keys():
            first_domain = '({},{})-1'.format(*configurations[0])
            domains = ['7', '8', first_domain, '11']
            simulation.strand_types['left_start'] = Strand(domains, False)
            simulation.strand_types['left_start_strip'] = Strand(domains, True)
    else:
        left_instruction = ['left_start_extended']
        if 'left_start_extended' not in simulation.strand_types.keys():
            domains = ['7', '8']
            for i in range(1, 3):
                domains.append('({},{})-{}'.format(*configurations[0], i))
            domains.append('11')
            simulation.strand_types['left_start_extended'] = Strand(domains, False)
            simulation.strand_types['left_start_extended_strip'] = Strand(domains, True)
    for i in range(0, index - 1):
        current_configuration = configurations[i]
        next_configuration = configurations[i + 1]
        left_strand_name = '({},{})_left'.format(*current_configuration)
        if left_strand_name not in simulation.strand_types.keys():
            left_strand_domains = ['({},{})-2'.format(*current_configuration),
                                   '({},{})-1'.format(*next_configuration), '11']
            color = simulation.strand_types['({},{})_full'.format(*current_configuration)].color
            simulation.strand_types[left_strand_name] = Strand(left_strand_domains, False, color)
        left_instruction.append(left_strand_name)
    if index > 0:
        current_configuration = configurations[index - 1]
        next_configuration = configurations[index]
        left_extended_name = '({},{})_left_extended'.format(*current_configuration)
        if left_extended_name not in simulation.strand_types.keys():
            left_extended_domains = ['({},{})-2'.format(*current_configuration)]
            for i in range(1, 3):
                left_extended_domains.append('({},{})-{}'.format(*next_configuration, i))
            left_extended_domains.append('11')
            color = simulation.strand_types['({},{})_full'.format(*current_configuration)].color
            simulation.strand_types[left_extended_name] = Strand(left_extended_domains, False, color)
            left_extended_strip_name = '({},{})_left_extended_strip'.format(*current_configuration)
            simulation.strand_types[left_extended_strip_name] = Strand(left_extended_domains, True, color)
        left_instruction.append(left_extended_name)
    simulation.instructions.append(left_instruction)

    if 'one_left_start' not in simulation.strand_types.keys():
        first_domain = '({},{})-1'.format(*configurations[0])
        color = simulation.strand_types['symbol_345678'].color
        domains = ['3', '4', '5', '6', '7', '8', first_domain, '13']
        simulation.strand_types['one_left_start'] = Strand(domains, False, color)
        simulation.strand_types['one_left_strip'] = Strand(domains, True, color)

    # Add instructions for when left cell is zero
    next_index = configurations.index((next_state, '0')) \
        if (next_state, '0') in configurations else -1
    if next_index == -1:
        if 'temp_zero' not in simulation.strand_types.keys():
            color = simulation.strand_types['symbol_45'].color
            domains = ['4', '5', '6', '13']
            simulation.strand_types['temp_zero'] = Strand(domains, False, color)
            simulation.strand_types['temp_zero_strip'] = Strand(domains, True, color)
        simulation.instructions.append(['temp_zero'])
        simulation.instructions.append(['one_left_start'])
        simulation.instructions.append(['temp_zero_strip'])
        simulation.instructions.append(['symbol_45', 'symbol_678'])
    else:
        zero_instruction = []
        generate_left_cascade_instruction_strands(next_index, num_configurations, configurations, zero_instruction,
                                                  simulation)
        if 'zero_symbol_displace' not in simulation.strand_types.keys():
            domains = ['2', '3', '4', '5', '6', '12']
            color = simulation.strand_types['symbol_45'].color
            simulation.strand_types['zero_symbol_displace'] = Strand(domains, False, color)
        zero_instruction.append('zero_symbol_displace')
        simulation.instructions.append(zero_instruction)
        simulation.instructions.append(['one_left_start'])  # in case left cell is 1
        generate_left_final_instruction_strands(next_index, num_configurations, configurations, simulation)

    # Add instructions for when left cell is one
    next_index = configurations.index((next_state, '1')) \
        if (next_state, '1') in configurations else -1
    simulation.instructions.append(['one_left_strip'])
    if next_index == -1:
        simulation.instructions.append(['symbol_345678'])
    else:
        if 'temp_one' not in simulation.strand_types.keys():
            color = simulation.strand_types['symbol_345678'].color
            first_domain = '({},{})-1'.format(*configurations[0])
            domains = ['4', '5', '6', '7', '8', first_domain, '13']
            simulation.strand_types['temp_one'] = Strand(domains, False, color)

        simulation.instructions.append(['temp_one'])

        one_instruction = []
        generate_left_cascade_instruction_strands(next_index, num_configurations, configurations, one_instruction,
                                                  simulation)
        if 'one_symbol_displace' not in simulation.strand_types.keys():
            domains = ['2', '3', '12']
            color = simulation.strand_types['symbol_12'].color
            simulation.strand_types['one_symbol_displace'] = Strand(domains, False, color)
        one_instruction.append('one_symbol_displace')
        simulation.instructions.append(one_instruction)
        generate_left_final_instruction_strands(next_index, num_configurations, configurations, simulation)

    # Add instructions for when left cell is blank
    next_index = configurations.index((next_state, blank_symbol)) \
        if (next_state, blank_symbol) in configurations else -1
    if index > 0:
        simulation.instructions.append(['left_start_strip'])
    else:
        simulation.instructions.append(['left_start_extended_strip'])

    if next_index == -1:
        simulation.instructions.append(['symbol_78'])
    else:
        simulation.instructions.append(['symbol_right'])
        simulation.instructions.append(['symbol_right_strip'])
        blank_instruction = []
        generate_left_cascade_instruction_strands(next_index, num_configurations, configurations, blank_instruction,
                                                  simulation)
        simulation.instructions.append(blank_instruction)
        generate_left_final_instruction_strands(next_index, num_configurations, configurations, simulation)

    # Replace current cell top strands
    if index > 0:
        simulation.instructions.append(['({},{})_left_extended_strip'.format(*configurations[index - 1])])
    else:
        simulation.instructions.append(['left_start_extended_strip'])

    replacement_instructions = []
    for i in range(0, index):
        replacement_instructions.append('({},{})_full'.format(*configurations[i]))
    generate_right_cascade_instruction_strands(index + 1, num_configurations, configurations, replacement_instructions,
                                               simulation)
    simulation.instructions.append(replacement_instructions)

    right_strip_instructions = []
    for i in range(index + 1, num_configurations):
        current_configuration = configurations[i]
        right_strip_instructions.append('({},{})_right_strip'.format(*current_configuration))
    right_strip_instructions.append('symbol_right_strip')
    simulation.instructions.append(right_strip_instructions)

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

    simulation.instructions.append(last_one_instructions)


def generate_left_cascade_instruction_strands(next_index, num_configurations, configurations, instruction_list,
                                              simulation):
    for i in range(next_index, num_configurations):
        current_configuration = configurations[i]
        next_configuration = configurations[i + 1] if i < num_configurations - 1 else -1
        left_strand_name = '({},{})_left'.format(*current_configuration)
        if left_strand_name not in simulation.strand_types.keys():
            left_strand_domains = ['({},{})-2'.format(*current_configuration)]
            if next_configuration == -1:
                left_strand_domains.append('1')
            else:
                left_strand_domains.append('({},{})-{}'.format(*next_configuration, 1))
            left_strand_domains.append('11')
            color = simulation.strand_types['({},{})_full'.format(*current_configuration)].color
            simulation.strand_types[left_strand_name] = Strand(left_strand_domains, False, color)
        instruction_list.append(left_strand_name)


def generate_left_final_instruction_strands(next_index, num_configurations, configurations, simulation):
    new_cell_instructions = ['({},{})_final'.format(*configurations[next_index])]
    for i in range(next_index + 1, num_configurations):
        new_cell_instructions.append('({},{})_full'.format(*configurations[i]))
    new_cell_instructions.append('symbol_covered')
    simulation.instructions.append(new_cell_instructions)


def create_right_instruction(index, num_configurations, transition, transition_data, blank_symbol, simulation):
    next_state = transition['R']
    write_symbol = transition['write']
    configurations = list(transition_data.keys())

    # Displace all strands until right end of cell
    right_instruction = []

    generate_right_cascade_instruction_strands(index + 1, num_configurations, configurations, right_instruction,
                                               simulation)
    simulation.instructions.append(right_instruction)

    # Strip out strands from previous instruction
    right_strip_instructions = []
    for i in range(index + 1, num_configurations):
        current_configuration = configurations[i]
        right_strip_instructions.append('({},{})_right_strip'.format(*current_configuration))
    right_strip_instructions.append('symbol_right_strip')

    simulation.instructions.append(right_strip_instructions)

    # Replace current cell top strands
    replacement_instructions = []
    for current_configuration in configurations:
        replacement_instructions.append('({},{})_full'.format(*current_configuration))

    if write_symbol == blank_symbol:
        replacement_instructions.append('symbol_123')
        if 'symbol_45(13)' not in simulation.strand_types.keys():
            domains = ['4', '5', '13']
            simulation.strand_types['symbol_45(13)'] = Strand(domains, False,
                                                              simulation.strand_types['symbol_456'].color)
            simulation.strand_types['symbol_45(13)_strip'] = Strand(domains, True,
                                                                    simulation.strand_types['symbol_456'].color)
            domains = ['6', '7', '13']
            simulation.strand_types['symbol_67(13)'] = Strand(domains, False,
                                                              simulation.strand_types['symbol_78'].color)
            simulation.strand_types['symbol_67(13)_strip'] = Strand(domains, True,
                                                                    simulation.strand_types['symbol_78'].color)
        replacement_instructions.extend(['symbol_45(13)', 'symbol_67(13)'])
    elif write_symbol == '0':
        replacement_instructions.append('symbol_123')
        replacement_instructions.append('symbol_45')
        if 'symbol_67' not in simulation.strand_types.keys():
            simulation.strand_types['symbol_67'] = Strand(['6', '7'], False,
                                                          simulation.strand_types['symbol_678'].color)
        replacement_instructions.append('symbol_67')
    else:
        replacement_instructions.append('symbol_12')
        if 'symbol_34567' not in simulation.strand_types.keys():
            simulation.strand_types['symbol_34567'] = Strand(['3', '4', '5', '6', '7'], False,
                                                             simulation.strand_types['symbol_345678'].color)
        replacement_instructions.append('symbol_34567')

    simulation.instructions.append(replacement_instructions)

    # Displace right cell contents
    next_right_instructions = []
    if '({},{})_next_right'.format(*configurations[0]) not in simulation.strand_types.keys():
        for i in range(len(configurations)):
            domains = ['8' if i == 0 else '({},{})-2'.format(*configurations[i - 1])]
            current_configuration = configurations[i]
            domains.append('({},{})-1'.format(*current_configuration))
            domains.append('11')
            strand_name = '({},{})_next_right'.format(*current_configuration)
            color = simulation.strand_types['({},{})_full'.format(*current_configuration)].color
            simulation.strand_types[strand_name] = Strand(domains, False, color)
            strand_strip_name = '({},{})_next_right_strip'.format(*current_configuration)
            simulation.strand_types[strand_strip_name] = Strand(domains, True, color)

        last_configuration = configurations[num_configurations - 1]
        domains = ['({},{})-2'.format(*last_configuration), '1', '2', '11']
        simulation.strand_types['next_right_1'] = Strand(domains, False)
        simulation.strand_types['next_right_1_strip'] = Strand(domains, True)
        simulation.strand_types['next_right_2'] = Strand(['3', '4', '5', '11'], False)
        simulation.strand_types['next_right_2_strip'] = Strand(simulation.strand_types['next_right_2'].domains, True)

    for current_configuration in configurations:
        next_right_instructions.append('({},{})_next_right'.format(*current_configuration))
    next_right_instructions.append('next_right_1')
    next_right_instructions.append('next_right_2')

    simulation.instructions.append(next_right_instructions)

    # Add instructions for when right cell is blank
    next_index = configurations.index((next_state, blank_symbol)) \
        if (next_state, blank_symbol) in configurations else -1
    first_configuration = configurations[0]
    if 'right_blank' not in simulation.strand_types.keys():
        domains = ['4', '5', '6', '12']
        simulation.strand_types['right_blank'] = Strand(domains, False)
        simulation.strand_types['right_blank_strip'] = Strand(domains, True)
        domains = ['8']
        for i in range(1, 3):
            domains.append('({},{})-{}'.format(*first_configuration, i))
        domains.append('12')
        color = simulation.strand_types['({},{})_full'.format(*first_configuration)].color
        simulation.strand_types['8_right_plug'] = Strand(domains, False, color)
        simulation.strand_types['8_right_plug_strip'] = Strand(domains, True, color)
    blank_instructions = ['8_right_plug']
    for i in range(1, len(configurations)):
        current_configuration = configurations[i]
        strand_type = '({},{})_full' if i != next_index else '({},{})_final'
        blank_instructions.append(strand_type.format(*current_configuration))
    blank_instructions.append('symbol_123')
    blank_instructions.append('right_blank')
    simulation.instructions.append(blank_instructions)
    simulation.instructions.append(['right_blank_strip'])
    simulation.instructions.append(['symbol_456' if next_index == -1 else 'symbol_covered'])
    simulation.instructions.append(['8_right_plug_strip'])
    last_blank_instruction = []
    if write_symbol == blank_symbol:
        last_blank_instruction.extend(['symbol_456', 'symbol_78'])
    elif write_symbol == '0':
        last_blank_instruction.append('symbol_678')
    else:
        last_blank_instruction.append('symbol_345678')
    strand_type = '({},{})_final' if next_index == 0 else '({},{})_full'
    last_blank_instruction.append(strand_type.format(*first_configuration))
    simulation.instructions.append(last_blank_instruction)

    # Add instructions for when right cell is zero
    next_index = configurations.index((next_state, '0')) \
        if (next_state, '0') in configurations else -1
    simulation.instructions.append(['next_right_2_strip'])
    if next_index == -1:
        simulation.instructions.append(['symbol_123', 'symbol_45'])
    else:
        simulation.instructions.append(['symbol_covered'])
    zero_instructions = ['8_right_plug']
    for i in range(1, len(configurations)):
        current_configuration = configurations[i]
        strand_type = '({},{})_full' if i != next_index else '({},{})_final'
        zero_instructions.append(strand_type.format(*current_configuration))
    simulation.instructions.append(zero_instructions)
    simulation.instructions.append(['8_right_plug_strip'])
    last_zero_instruction = []
    if write_symbol == blank_symbol:
        last_zero_instruction.extend(['symbol_456', 'symbol_78'])
    elif write_symbol == '0':
        last_zero_instruction.append('symbol_678')
    else:
        last_zero_instruction.append('symbol_345678')
    strand_type = '({},{})_final' if next_index == 0 else '({},{})_full'
    last_zero_instruction.append(strand_type.format(*first_configuration))
    simulation.instructions.append(last_zero_instruction)

    # Add instructions for when right cell is one
    next_index = configurations.index((next_state, '1')) \
        if (next_state, '1') in configurations else -1
    one_instructions = []
    for current_configuration in configurations:
        one_instructions.append('({},{})_next_right_strip'.format(*current_configuration))
    one_instructions.append('next_right_1_strip')
    simulation.instructions.append(one_instructions)
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
    simulation.instructions.append(next_one_instruction)

    if next_index == -1:
        simulation.instructions.append(['symbol_12', 'symbol_345678'])
    else:
        simulation.instructions.append(['symbol_covered'])


def generate_right_cascade_instruction_strands(index, num_configurations, configurations, instruction_list, simulation):
    for i in range(index, num_configurations):
        previous_configuration = configurations[i - 1] if i > 0 else -1
        current_configuration = configurations[i]
        right_strand_name = '({},{})_right'.format(*current_configuration)
        if right_strand_name not in simulation.strand_types.keys():
            if previous_configuration == -1:
                right_strand_domains = ['8']
            else:
                right_strand_domains = ['({},{})-2'.format(*previous_configuration)]
            right_strand_domains.append('({},{})-1'.format(*current_configuration))
            right_strand_domains.append('11')
            color = simulation.strand_types['({},{})_full'.format(*current_configuration)].color
            simulation.strand_types[right_strand_name] = Strand(right_strand_domains, False, color)
            right_strip_name = '({},{})_right_strip'.format(*current_configuration)
            simulation.strand_types[right_strip_name] = Strand(right_strand_domains, True, color)
        instruction_list.append(right_strand_name)

    if 'symbol_right' not in simulation.strand_types.keys():
        last_configuration = configurations[num_configurations - 1]
        right_strand_domains = ['({},{})-2'.format(*last_configuration)]
        for i in range(1, 8):
            right_strand_domains.append(str(i))
        right_strand_domains.append('11')
        color = simulation.strand_types['symbol_covered'].color
        simulation.strand_types['symbol_right'] = Strand(right_strand_domains, False, color)
        simulation.strand_types['symbol_right_strip'] = Strand(right_strand_domains, True, color)

    instruction_list.append('symbol_right')
