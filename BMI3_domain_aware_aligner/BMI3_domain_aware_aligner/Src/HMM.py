from Src.multi_alignment import *
import numpy as np


def simulate_data(identifier, number, alphabet, length, mutation, deletion):
    '''
    Generate simulate data using given patern
    '''
    import random
    aalist = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
              'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    def random_seq_generator(size, word_list):
        res = ''.join(random.choice(word_list) for _ in range(size))
        return res

    def randomly_generate_mutation(seq, mutation):
        length = len(seq)
        number_of_mutation = int(np.random.normal(
            loc=mutation[0], scale=mutation[1]) * length)
        points = [random.choice(list(range(0, length)))
                  for _ in range(number_of_mutation)]
        seql = list(seq)
        for point in points:
            seql[point] = random.choice(aalist)
        return ''.join(seql)

    def randomly_generate_deletion(seq, deletion):
        length = len(seq)
        number_of_deletion = int(np.random.normal(
            loc=deletion[0], scale=deletion[1]) * length)
        points = [random.choice(list(range(0, length)))
                  for _ in range(number_of_deletion)]
        seql = list(seq)
        for point in points:
            seql[point] = '-'
        return ''.join(seql).replace('-', '')

    res_dict = {'domains': {}, 'sequences': {}}
    domains = {}
    for i in range(len(alphabet)):
        _name = alphabet[i]
        _length = length[i]
        seq = random_seq_generator(_length, aalist)
        mutated_seq = randomly_generate_mutation(seq, mutation)
        domains[_name] = seq

    sequences = {}
    seq_info = {}
    for i in range(len(identifier)):
        id = identifier[i]
        number_i = number[i]
        sequences[id] = []
        seq_info[id] = []
        for n in range(number_i):
            # generate n sequences of identifier id
            sections = []
            current_index = 0
            section_dict = {}
            counter = {}
            for s in list(id):
                if s in counter:
                    counter[s] += 1
                else:
                    counter[s] = 0
                seq = domains[s]
                mutated_seq = randomly_generate_mutation(seq, mutation)
                mutated_seq = randomly_generate_deletion(mutated_seq, deletion)
                sections.append(mutated_seq)
                if not s in section_dict:
                    section_dict[s] = [current_index,
                                       current_index+len(mutated_seq)-1]
                else:
                    section_dict[s+str(counter[s])] = [current_index,
                                                       current_index+len(mutated_seq)-1]
                current_index += len(mutated_seq)

            seq_info[id].append(section_dict)
            sequences[id].append(''.join(sections))

    def rmlinker(seq_info, linkers):
        for i in seq_info:
            for j in range(len(seq_info[i])):
                for k in linkers:
                    if k in seq_info[i][j]:
                        del seq_info[i][j][k]
        return seq_info

    seq_info = rmlinker(seq_info, ['q', 'w', 'e', 'r', 't'])

    return domains, sequences, seq_info


def run_simulated_data():
    example_pattern = {
        # 'identifier'   :['qAwBeCr','qAwCeBr','qAwCeDr','wBeCr'],
        'identifier': ['qAwCeBr', 'qAwCeDr', 'wCeBr'],
        'number': [11, 6, 3, 7],
        'alphabet': ['A', 'B', 'C', 'D', 'q', 'w', 'e', 'r'],
        'length': [30, 40, 25, 45, 7, 9, 5, 7],
        'mutation': (0.01, 1),
        'deletion': (0.02, 0.05),
    }
    def simulate_data(identifier, number, alphabet, length, mutation, deletion):
        '''
        Generate simulate data using given patern
        '''
        import random
        aalist = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

        def random_seq_generator(size, word_list):
            res = ''.join(random.choice(word_list) for _ in range(size))
            return res

        def randomly_generate_mutation(seq, mutation):
            length = len(seq)
            number_of_mutation = int(np.random.normal(
                loc=mutation[0], scale=mutation[1]) * length)
            points = [random.choice(list(range(0, length))) for _ in range(number_of_mutation)]
            seql = list(seq)
            for point in points:
                seql[point] = random.choice(aalist)
            return ''.join(seql)

        def randomly_generate_deletion(seq, deletion):
            length = len(seq)
            number_of_deletion = int(np.random.normal(
                loc=deletion[0], scale=deletion[1]) * length)
            points = [random.choice(list(range(0, length))) for _ in range(number_of_deletion)]
            seql = list(seq)
            for point in points:
                seql[point] = '-'
            return ''.join(seql).replace('-', '')

        res_dict = {'domains': {}, 'sequences': {}}
        domains = {}
        for i in range(len(alphabet)):
            _name = alphabet[i]
            _length = length[i]
            seq = random_seq_generator(_length, aalist)
            mutated_seq = randomly_generate_mutation(seq, mutation)
            domains[_name] = seq

        sequences = {}
        seq_info = {}
        for i in range(len(identifier)):
            id = identifier[i]
            number_i = number[i]
            sequences[id] = []
            seq_info[id] = []
            for n in range(number_i):
                # generate n sequences of identifier id
                sections = []
                current_index = 0
                section_dict = {}
                for s in list(id):
                    seq = domains[s]
                    mutated_seq = randomly_generate_mutation(seq, mutation)
                    mutated_seq = randomly_generate_deletion(mutated_seq, deletion)
                    sections.append(mutated_seq)
                    section_dict[s] = [current_index, current_index + len(mutated_seq) - 1]
                    current_index += len(mutated_seq)
                seq_info[id].append(section_dict)
                sequences[id].append(''.join(sections))

        return domains, sequences, seq_info
    domains, sequences, seq_info = simulate_data(**example_pattern)
    def rmlinker(seq_info, linkers):
        for i in seq_info:
            for j in range(len(seq_info[i])):
                for k in linkers:
                    if k in seq_info[i][j]:
                        del seq_info[i][j][k]
    rmlinker(seq_info, ['q', 'w', 'e', 'r'])
    #######################################

    input_data = {}
    all_domains = {}
    for key in sequences:
        for pos in range(len(sequences[key])):
            input_data[key+str(pos)] = sequences[key][pos]
            all_domains[key+str(pos)] = seq_info[key][pos]
    new_all_domains = {}
    for key in all_domains:
        curr_arr = []
        dictionary = all_domains[key]
        for domain_key in dictionary:
            domain_start_end = dictionary[domain_key]
            curr_arr.append([domain_key, domain_start_end[0], domain_start_end[1]])
        new_all_domains[key] = curr_arr
    domain_alignment_result = domain_aware_greedy_MSA(new_all_domains, input_data)
    for structure in domain_alignment_result:
        print('-------------------------------------------------------')
        print(
            'Alignment result for ' + str(len(domain_alignment_result[structure]['seqs'])) + ' sequences of ' + structure +
            ' structure: ')
        for alignment in domain_alignment_result[structure]['category_alignment']:
            print(alignment)
    return domain_alignment_result, new_all_domains

    
def transform_to_standard_input_format(domains, sequences, seq_info):
    '''
    Transform simulated data into the format that can be directory used in subsequent analysis
    '''
    input_data = {}
    all_domains = {}
    for key in sequences:
        for pos in range(len(sequences[key])):
            input_data[key+str(pos)] = sequences[key][pos]
            all_domains[key+str(pos)] = seq_info[key][pos]
    new_all_domains = {}
    for key in all_domains:
        curr_arr = []
        dictionary = all_domains[key]
        for domain_key in dictionary:
            domain_start_end = dictionary[domain_key]
            curr_arr.append(
                [domain_key, domain_start_end[0], domain_start_end[1]])
        new_all_domains[key] = curr_arr
        # input_data, new_all_domains are used in func:
    return input_data, new_all_domains






def construct_HMM(theta, Alphabet, Alignment):
    '''
    Construct a HMM profile using aligned sequence.
    :param theta: Threshold of deciding whether a position should be a insertion or deletion.
    :param Alphabet: List of domain symbols appear in alignment.
    :param Alignment: Aligned sequences of domains. List of str format, strs should be of the same length.
                      
    :return StateIndices: List of tuple represent states of HMM in order of time & 'I,M,D'.
    :return Transition: Dictionary of transition rate between two state.
    :return Emission: Dictionary of emission rate of states. For D state, emission_rate=1.
    '''

    ### Defining tool functions ###

    def count_one_alignment(m, j, k):
        '''
        Used to count alphabetical characters in specified column of alignment
        m: Number of sequences
        j: The column to be counted
        k: Number of symbols in Alphabet
        return a nested list whose row represent position in alignment, column represent order in Alphabet
        '''
        Counts = [0]*k
        for i in range(m):
            if Alignment[i][j] in Alphabet:
                Counts[Alphabet.index(Alignment[i][j])] += 1
        return Counts

    def create_state_indices(valid_position):
        '''
        :param valid_position list showing position conservation
        :Given a list showing position is insertion or not, this function return a list of states in HMM
        '''
        index = 0
        state_indices = [('S', index), ('I', index)]
        for p in valid_position:
            if not p:
                continue
            index += 1
            state_indices.append(('M', index))
            state_indices.append(('D', index))
            state_indices.append(('I', index))
        index += 1
        state_indices.append(('E', index))
        return state_indices

    def create_state_counts(StateIndices):
        '''
        This function create a frame work illustrating state-state transfer.
        '''
        state_counts = {}
        # Get every pair of states, examine whether there is a one-way transition relationship between them.
        # (s,i) is always regard as the frontal one and (ns,ni) is regard as the posterior one.
        for s, i in StateIndices:
            for ns, ni in StateIndices:
                if s == 'S':
                    if ns == 'I' and ni == 0:
                        state_counts[(s, i), (ns, ni)] = 0
                    elif ns in ['M', 'D'] and ni == 1:
                        state_counts[(s, i), (ns, ni)] = 0
                elif s == 'I':
                    # Insertion state is special as it can produce self transition, so the time of HMM chain
                    # wouldn't change, so i should be equal to ni when self transition happen. But for transition
                    # of changing to Deletion and Match state, ni should be equal to i+1.
                    if i == 0:
                        if ns == 'I' and ni == 0:
                            state_counts[(s, i), (ns, ni)] = 0
                        elif ns in ['M', 'D'] and ni == 1:
                            state_counts[(s, i), (ns, ni)] = 0
                    else:
                        if (ns == 'I' and ni == i) or (ns in ['M', 'D', 'E'] and ni == i+1):
                            state_counts[(s, i), (ns, ni)] = 0
                elif s in ['M', 'D']:
                    # If current state is one of Matches or Deltion and its going to transfer to state Insertion,
                    # the time of HMM chain wouln't change either.
                    if (ns == 'I' and ni == i) or (ns in ['M', 'D', 'E'] and ni == i+1):
                        state_counts[(s, i), (ns, ni)] = 0
                else:
                    assert(
                        s == 'E'), 'An unexception error occur, please check your input!'
        return state_counts

    def create_emission_counts(StateIndices, Alphabet):
        '''
        This function create a frame work illustrating emission rate of states.
        '''
        emission_counts = {}
        for state in StateIndices:
            for symbol in Alphabet:
                emission_counts[(state, symbol)] = 0
        return emission_counts

    def create_state_frequencies(StateCounts, StateIndices):
        for i in StateCounts.keys():
            StateCounts[i] += m/10  # use global varible m
        Totals = {i: 0 for i in StateIndices}

        for key, count in StateCounts.items():
            Totals[key[0]] += count

        result = {}
        for key, count in StateCounts.items():
            result[key] = count/Totals[key[0]] if Totals[key[0]] > 0 else 0

        return result

    def create_emission_frequencies(EmissionCounts, StateIndices):
        for i in EmissionCounts.keys():
            if i[0][0] in ['M', 'I']:
                EmissionCounts[i] += k/100  # use global varible k
        Totals = {i: 0 for i in StateIndices}

        for key, count in EmissionCounts.items():
            Totals[key[0]] += count

        result = {}
        for key, count in EmissionCounts.items():
            result[key] = count/Totals[key[0]] if Totals[key[0]] > 0 else 0

        return result

    ### Main body of building HMM profile function ###

    # Define commonly used variables
    global k, m, n
    k, m, n = len(Alphabet), len(Alignment), len(Alignment[0])
    # make sure that all sequences' lenght is the same
    for s in Alignment:
        assert(n == len(s)), 'All sequece input should be of the same length!'

    # count number of symbols in each column
    counts = [count_one_alignment(m, j, k) for j in range(n)]
    # Detect whether or not 'insertion_rate' in positions are larger than the given threshold
    valid_position = [sum(c) > (1-theta)*m for c in counts]

    # Create state indices.
    # Initialize framework of transition rate and emission rate
    Inserts, Matches, Deletes = [[[] for _ in range(sum(valid_position)+1)]]*3
    StateIndices = create_state_indices(valid_position)
    StateCounts = create_state_counts(StateIndices)
    EmissionCounts = create_emission_counts(StateIndices, Alphabet)

    # Assign count to each state in HMM and also with the emission infomation
    for seq in Alignment:
        prev_state_type = 'S'
        state_type_history = [prev_state_type]
        time = 0
        # First, we need to verify how the n symbols in each aligned sequence pass through HMM.
        # Variable state_type_history is used for recording them.
        for i in range(n):
            symbol = seq[i]
            if valid_position[i]:
                time += 1
                if symbol in Alphabet:
                    # Matches[time].append((prev_state_type, symbol))
                    prev_state_type = 'M'
                    state_type_history.append(prev_state_type)
                elif symbol == '-':
                    # Deletes[time].append((prev_state_type, symbol))
                    prev_state_type = 'D'
                    state_type_history.append(prev_state_type)
            else:
                if symbol in Alphabet:
                    # Inserts[time].append((prev_state_type,symbol))
                    prev_state_type = 'I'
                    state_type_history.append(prev_state_type)
                elif symbol == '-':
                    pass
        state_type_history.append('E')
        print(f'Counting {seq} {"".join(state_type_history)}')

        # Next, we fill each aligned sequence into StateCounts and EmissionCounts dictionary,
        # which will be sebsequently used to calculate transition rate and emission rate.
        prev_state = None
        time, seq_index = 0, 0
        for state in state_type_history:
            while seq_index < len(seq)-1 and seq[seq_index] == '-':
                seq_index += 1
            if state in ['M', 'D', 'E']:
                time += 1

            if prev_state:
                StateCounts[prev_state, (state, time)] += 1
                if state in ['M', 'I']:
                    EmissionCounts[(state, time), seq[seq_index]] += 1
                    seq_index += 1
            prev_state = (state, time)

    Transition = create_state_frequencies(StateCounts, StateIndices)
    Emission = create_emission_frequencies(EmissionCounts, StateIndices)

    return StateIndices, Transition, Emission


def align_to_HMM(xs, Alphabet, States, Transition, Emission):
    '''
     start state and end state is ('S',0) & ('E',times) by default
    :param xs: sequence of domains
    :param Alphabet: dictionary all domain appear in sequence
    :param States: dictionary of all state of HMM
    :param Transition: dictionary of transition rate between two state
    :param Emission: dictionary of emission rate of states
    
    :return: viterbi alignment result and its most probable hidden path in HMM
    '''

    # Define commonly used variables
    length = len(xs)
    # time represent the section or column of hidden markov chain
    times = States[-1][1]-1

    # Initialize prev_states dictionary
    prev_states = {}
    for state in States:
        # state in format: < ('X',digit) >,
        # where 'X' represent state's typr and digit represent timepoint
        state_type, time = state
        if state_type in ['D', 'M', 'E']:
            # As deletion, match and end state will bring a +1 increase of time along
            # the hidden markov chain, their previous state time is their time-1.
            prev_states[state] = [('D', time-1), ('I', time-1), ('M', time-1)]
        elif state_type == 'I':
            # While insertion state just add a 'gap' into the modol, it will not increase time of state.
            prev_states[state] = [('D', time), ('I', time), ('M', time)]

    # These three state need to be reset manually as they are connect to S0 state but not infered regular states.
    prev_states[('I', 0)] = [('I', 0)]
    prev_states[('M', 1)] = [('I', 0)]
    prev_states[('D', 1)] = [('I', 0)]

    # Initialize probability matrix for dynamic programming
    dp, memory = {}, {}

    for state in States[:-1]:
        dp[state] = [0.0]*(length+1)
        memory[state] = [[] for _ in range(length+1)]

    dp[('S', 0)][0] = 1.0
    dp[('I', 0)][1] = Transition[(('S', 0), ('I', 0))]*1
    memory[('I', 0)][1] = (('S', 0), xs[0])
    dp[('M', 1)][1] = Transition[(('S', 0), ('M', 1))] * \
        Emission[(('M', 1), xs[0])]*1
    memory[('M', 1)][1] = (('S', 0), xs[0])

    for t in range(times):  # first column only contain deletion state
        dp[('D', 1)][0] = Transition[(('S', 0), ('D', 1))]
        for i in range(2, times+1):
            dp[('D', i)][0] = Transition[(('D', i-1), ('D', i))]

    # Start calculation best situation of every state just considering its previous states
    def get_best_prev(l, state, dp, memory, Transition, Emission):
        current_s = xs[l-1]
        state_type, time = state
        candidates = []
        if state_type == 'M':
            for prev_state in prev_states[state]:
                transition_rate = Transition[(prev_state, state)]
                candidates.append(
                    (prev_state, dp[prev_state][l-1]*transition_rate*Emission[(state, current_s)]))

        elif state_type == 'I':
            for prev_state in prev_states[state]:
                transition_rate = Transition[(prev_state, state)]
                candidates.append(
                    (prev_state, dp[prev_state][l-1]*transition_rate*Emission[(state, current_s)]))
        elif state_type == 'D':
            for prev_state in prev_states[state]:
                transition_rate = Transition[(prev_state, state)]
                candidates.append(
                    (prev_state, dp[prev_state][l]*transition_rate))

        best_prev, dp_score = max(candidates, key=lambda x: x[1])
        dp[state][l] = dp_score

        emit = '-' if state_type == 'D' else current_s
        memory[state][l] = (best_prev, emit)
        return dp, memory

    for l in range(1, length+1):
        for state in States[1:-1]:
            if l == 1 and (state in [('I', 0), ('M', 1)]):
                continue
            dp, memory = get_best_prev(
                l, state, dp, memory, Transition, Emission)

    # for end state
    end_state = States[-1]
    candidates = []
    for prev_state in States[-4:-1]:
        candidates.append(
            (prev_state, dp[prev_state][length]*Transition[(prev_state, end_state)]))

    last, dp_score = max(candidates, key=lambda x: x[1])
    l = length
    path = [(end_state, l+1)]
    res = ''
    while memory[last][l]:
        path.append((last, l))
        last, emit = memory[last][l]
        res = emit+res
        if last[0] != 'D':
            l = l-1

    path.append((('S', 0), 0))
    path.reverse()
    return dp, memory, res, path


def construct_SWAP_HMM(theta, Alphabet, Alignment):
    '''
    Construct a HMM profile using aligned sequence.
    :param theta: Threshold of deciding whether a position should be a insertion or deletion.
    :param Alphabet: List of domain symbols appear in alignment.
    :param Alignment: Aligned sequences of domains. List of str format, strs should be of the same length.
                      
    :return StateIndices: List of tuple represent states of HMM in order of time & 'I,M,D'.
    :return Transition: Dictionary of transition rate between two state.
    :return Emission: Dictionary of emission rate of states. For D state, emission_rate=1.
    '''

    ### Defining tool functions ###

    def count_one_alignment(m, j, k):
        '''
        Used to count alphabetical characters in specified column of alignment
        m: Number of sequences
        j: The column to be counted
        k: Number of symbols in Alphabet
        return a nested list whose row represent position in alignment, column represent order in Alphabet
        '''
        Counts = [0]*k
        for i in range(m):
            if Alignment[i][j] in Alphabet:
                Counts[Alphabet.index(Alignment[i][j])] += 1
        return Counts

    def create_state_indices(valid_position):
        '''
        :param valid_position list showing position conservation
        :Given a list showing position is insertion or not, this function return a list of states in HMM
        '''
        index = 0
        state_indices = [('S', index), ('I', index)]
        for p in valid_position:
            if not p:
                continue
            index += 1
            state_indices.append(('M', index))
            state_indices.append(('D', index))
            # (need comment)
            if index >= 2:
                state_indices.append(('SW', index))
            state_indices.append(('I', index))
        index += 1
        state_indices.append(('E', index))
        return state_indices

    def create_state_counts(StateIndices):
        '''
        This function create a frame work illustrating state-state transfer.
        '''
        state_counts = {}
        # Get every pair of states, examine whether there is a one-way transition relationship between them.
        # (s,i) is always regard as the frontal one and (ns,ni) is regard as the posterior one.
        for s, i in StateIndices:
            for ns, ni in StateIndices:
                if s == 'S':
                    if ns == 'I' and ni == 0:
                        state_counts[(s, i), (ns, ni)] = 0
                    elif ns in ['M', 'D'] and ni == 1:
                        state_counts[(s, i), (ns, ni)] = 0
                    # (need comment)
                    elif ns == 'SW' and ni == 2:
                        state_counts[(s, i), (ns, ni)] = 0
                elif s == 'I':
                    # Insertion state is special as it can produce self transition, so the time of HMM chain
                    # wouldn't change, so i should be equal to ni when self transition happen. But for transition
                    # of changing to Deletion and Match state, ni should be equal to i+1.
                    if i == 0:
                        if ns == 'I' and ni == 0:
                            state_counts[(s, i), (ns, ni)] = 0
                        elif ns in ['M', 'D'] and ni == 1:
                            state_counts[(s, i), (ns, ni)] = 0
                        # (need comment)
                        elif ns == 'SW' and ni == 2:
                            state_counts[(s, i), (ns, ni)] = 0
                    else:
                        if (ns == 'I' and ni == i) or (ns in ['M', 'D', 'E'] and ni == i+1):
                            state_counts[(s, i), (ns, ni)] = 0
                        # (need comment)
                        elif ns == 'SW' and ni == i+2:
                            state_counts[(s, i), (ns, ni)] = 0
                elif s in ['M', 'D']:
                    # If current state is one of Matches or Deltion and its going to transfer to state Insertion,
                    # the time of HMM chain wouln't change either.
                    if (ns == 'I' and ni == i) or (ns in ['M', 'D', 'E'] and ni == i+1):
                        state_counts[(s, i), (ns, ni)] = 0
                    # (need comment)*
                    elif ns == 'SW' and ni == i+2:
                        state_counts[(s, i), (ns, ni)] = 0
                elif s == 'SW':
                    # If current state is one of Matches or Deltion and its going to transfer to state Insertion,
                    # the time of HMM chain wouln't change either.
                    if (ns == 'I' and ni == i) or (ns in ['M', 'D', 'E'] and ni == i+1):
                        state_counts[(s, i), (ns, ni)] = 0
                    # (need comment)*
                    elif ns == 'SW' and ni == i+2:
                        state_counts[(s, i), (ns, ni)] = 0
                else:
                    assert(
                        s == 'E'), 'An unexception error occur, please check your input!'
        return state_counts

    def create_emission_counts(StateIndices, Alphabet):
        '''
        This function create a frame work illustrating emission rate of states.
        '''
        emission_counts = {}
        for state in StateIndices:
            for symbol in Alphabet:
                emission_counts[(state, symbol)] = 0
        return emission_counts

    def create_state_frequencies(StateCounts, StateIndices):
        for i in StateCounts.keys():
            StateCounts[i] += m/10  # use global varible m
            # (need comment)**  explain why the pesudoCount is important for this model and why state SW was specially treated
            if i[1][0] == 'SW':
                StateCounts[i] += m/10
            if i[1][0] == 'I' and i[0][0] == 'I':
                StateCounts[i] -= m/15
        Totals = {i: 0 for i in StateIndices}

        for key, count in StateCounts.items():
            Totals[key[0]] += count

        result = {}
        for key, count in StateCounts.items():
            result[key] = count/Totals[key[0]] if Totals[key[0]] > 0 else 0

        return result

    def create_emission_frequencies(EmissionCounts, StateIndices):
        for i in EmissionCounts.keys():
            if i[0][0] in ['M', 'I']:
                EmissionCounts[i] += k/1000  # use global varible k
        Totals = {i: 0 for i in StateIndices}

        for key, count in EmissionCounts.items():
            Totals[key[0]] += count

        result = {}
        for key, count in EmissionCounts.items():
            result[key] = count/Totals[key[0]] if Totals[key[0]] > 0 else 0

        return result

    ### Main body of building HMM profile function ###

    # Define commonly used variables
    global k, m, n
    k, m, n = len(Alphabet), len(Alignment), len(Alignment[0])
    # make sure that all sequences' lenght is the same
    for s in Alignment:
        assert(n == len(s)), 'All sequece input should be of the same length!'

    # count number of symbols in each column
    counts = [count_one_alignment(m, j, k) for j in range(n)]
    # Detect whether or not 'insertion_rate' in positions are larger than the given threshold
    valid_position = [sum(c) > (1-theta)*m for c in counts]

    # Create state indices.
    # Initialize framework of transition rate and emission rate
    StateIndices = create_state_indices(valid_position)
    StateCounts = create_state_counts(StateIndices)
    EmissionCounts = create_emission_counts(StateIndices, Alphabet)

    # Assign count to each state in HMM and also with the emission infomation
    # In the profile construction step, we will not consider the swap state at first, but just initialize them
    for seq in Alignment:
        prev_state_type = 'S'
        state_type_history = [prev_state_type]
        time = 0
        # First, we need to verify how the n symbols in each aligned sequence pass through HMM.
        # Variable state_type_history is used for recording them.
        for i in range(n):
            symbol = seq[i]
            if valid_position[i]:
                time += 1
                if symbol in Alphabet:
                    # Matches[time].append((prev_state_type, symbol))
                    prev_state_type = 'M'
                    state_type_history.append(prev_state_type)
                elif symbol == '-':
                    # Deletes[time].append((prev_state_type, symbol))
                    prev_state_type = 'D'
                    state_type_history.append(prev_state_type)
            else:
                if symbol in Alphabet:
                    # Inserts[time].append((prev_state_type,symbol))
                    prev_state_type = 'I'
                    state_type_history.append(prev_state_type)
                elif symbol == '-':
                    pass
        state_type_history.append('E')
        # print(f'Counting {seq} {"".join(state_type_history)}')

        # Next, we fill each aligned sequence into StateCounts and EmissionCounts dictionary,
        # which will be sebsequently used to calculate transition rate and emission rate.
        prev_state = None
        time, seq_index = 0, 0
        for state in state_type_history:
            while seq_index < len(seq)-1 and seq[seq_index] == '-':
                seq_index += 1
            if state in ['M', 'D', 'E']:
                time += 1

            if prev_state:
                StateCounts[prev_state, (state, time)] += 1
                if state in ['M', 'I']:
                    EmissionCounts[(state, time), seq[seq_index]] += 1
                    seq_index += 1
            prev_state = (state, time)

    Transition = create_state_frequencies(StateCounts, StateIndices)
    Emission = create_emission_frequencies(EmissionCounts, StateIndices)

    return StateIndices, Transition, Emission


def align_to_SWAP_HMM(xs, Alphabet, States, Transition, Emission):
    '''
    (need comment)**
    :breief: start state and end state is ('S',0) & ('E',times) by default
    :param xs: sequence of domains
    :param Alphabet: dictionary all domain appear in sequence
    :param States: dictionary of all state of HMM
    :param Transition: dictionary of transition rate between two state
    :param Emission: dictionary of emission rate of states
    
    :return: viterbi alignment result and its most probable hidden path in HMM
    '''

    # Define commonly used variables
    length = len(xs)
    # time represent the section or column of hidden markov chain
    times = States[-1][1]-1

    # (need comment)
    def get_all_Matches(times, Alphabet, Emission):
        res = ''
        for i in range(1, times):
            tmp = []
            for s in Alphabet:
                tmp.append((Emission[(('M', i), s)], s))  # type: ignore
            m = max(tmp, key=lambda x: x[0])[1]
            res += m
        return res

    matchStates = get_all_Matches(times+1, Alphabet, Emission)

    # Initialize prev_states dictionary
    prev_states = {}
    for state in States:
        # state in format: < ('X',digit) >,
        # where 'X' represent state's typr and digit represent timepoint
        state_type, time = state
        if state_type in ['D', 'M', 'E']:
            # As deletion, match and end state will bring a +1 increase of time along
            # the hidden markov chain, their previous state time is their time-1.
            prev_states[state] = [('D', time-1), ('I', time-1), ('M', time-1)]
            if time >= 3:
                prev_states[state].append(('SW', time-1))
        elif state_type == 'I':
            # While insertion state just add a 'gap' into the modol, it will not increase time of state.
            prev_states[state] = [('D', time), ('I', time), ('M', time)]
            if time >= 2:
                prev_states[state].append(('SW', time))
        elif state_type == 'SW':
            # (need comment)*
            prev_states[state] = [('D', time-2), ('I', time-2), ('M', time-2)]
            if time >= 4:
                prev_states[state].append(('SW', time-2))

    # These four state need to be reset manually as they are connect to S0 state but not infered regular states.
    prev_states[('I', 0)] = [('I', 0)]
    prev_states[('M', 1)] = [('I', 0)]
    prev_states[('D', 1)] = [('I', 0)]
    prev_states[('SW', 2)] = [('I', 0)]

    # Initialize probability matrix for dynamic programming
    dp, memory = {}, {}

    for state in States[:-1]:
        dp[state] = [0.0]*(length+1)
        memory[state] = [[] for _ in range(length+1)]

    dp[('S', 0)][0] = 1.0
    # (need comment)*
    dp[('I', 0)][1] = Transition[(('S', 0), ('I', 0))]
    memory[('I', 0)][1] = (('S', 0), xs[0])
    dp[('M', 1)][1] = Transition[(('S', 0), ('M', 1))] * \
        Emission[(('M', 1), xs[0])]*1
    memory[('M', 1)][1] = (('S', 0), xs[0])
    # (need comment)

    if xs[:2][::-1] == matchStates[:2]:  # xs[l-2:l][::-1] == matchStates[time-2:time]
        memory[('SW', 2)][2] = (('S', 0), xs[0:2][::-1])
        dp[('SW', 2)][2] = Transition[(('S', 0), ('SW', 2))]

     # first column only contain deletion state
    dp[('D', 1)][0] = Transition[(('S', 0), ('D', 1))]
    memory[('D', 1)][0] = (('S', 0), '-')
    for i in range(2, times+1):
        dp[('D', i)][0] = Transition[(('D', i-1), ('D', i))]
        memory[('D', i)][0] = (('D', i-1), '-')
    # Start calculation best situation of every state just considering its previous states

    def get_best_prev(l, state, dp, memory, Transition, Emission):
        '''
        (need comment)*
        '''
        current_s = xs[l-1]
        state_type, time = state
        candidates = []
        '''
        # if state_type == 'M':
        #     for prev_state in prev_states[state]:
        #         transition_rate = Transition[(prev_state, state)]
        #         candidates.append(
        #             (prev_state, dp[prev_state][l-1]*transition_rate*Emission[(state, current_s)]))

        # elif state_type == 'I':
        #     for prev_state in prev_states[state]:
        #         transition_rate = Transition[(prev_state, state)]
        #         candidates.append(
        #             (prev_state, dp[prev_state][l-1]*transition_rate*Emission[(state, current_s)]))

        # elif state_type == 'D':
        #     for prev_state in prev_states[state]:
        #         transition_rate = Transition[(prev_state, state)]
        #         candidates.append(
        #             (prev_state, dp[prev_state][l]*transition_rate))

        # # (need comment)**
        # elif l >= 2 and state_type == 'SW':
        #     current_s = xs[l-2:l][::-1]
        #     for prev_state in prev_states[state]:
        #         transition_rate = Transition[(prev_state, state)]
        #         print(current_s)
        #         print(matchStates[time-2:time])
        #         if current_s == matchStates[time-2:time]:
        #             candidates.append((prev_state, dp[prev_state][l-2]*transition_rate))
        '''
        for prev_state in prev_states[state]:
            # if prev_state[0] == 'SW' and xs[l-2:l][::-1] != matchStates[time-2:time]:
            #     continue
            transition_rate = Transition[(prev_state, state)]
            if state_type == 'M':
                candidates.append(
                    (prev_state, dp[prev_state][l-1]*transition_rate*Emission[(state, current_s)]**2))
            elif state_type == 'I':
                candidates.append(
                    (prev_state, dp[prev_state][l-1]*transition_rate*Emission[(state, current_s)]))
            elif state_type == 'D':
                candidates.append(
                    (prev_state, dp[prev_state][l]*transition_rate))
            # (need comment)**
            elif l >= 2 and state_type == 'SW':
                current_s = xs[l-2:l][::-1]
                if current_s == matchStates[time-2:time]:
                    candidates.append(
                        (prev_state, dp[prev_state][l-2]*transition_rate))

        # (need comment)
        if candidates:
            best_prev, dp_score = max(candidates, key=lambda x: x[1])
            if dp_score > dp[state][l]:
                dp[state][l] = dp_score
                emit = '-' if state_type == 'D' else current_s
                memory[state][l] = (best_prev, emit)

        return dp, memory

    for l in range(1, length+1):
        for state in States[1:-1]:
            # (need comment)
            if (l == 1 and (state in [('I', 0), ('M', 1)])) or (l == 2 and state == 'SW'):
                continue
            dp, memory = get_best_prev(
                l, state, dp, memory, Transition, Emission)

    # for end state
    end_state = States[-1]
    candidates = []
    # (need comment)
    for prev_state in States[-5:-1]:
        candidates.append(
            (prev_state, dp[prev_state][length]*Transition[(prev_state, end_state)]))

    # print(candidates)
    last, dp_score = max(candidates, key=lambda x: x[1])
    l = length
    path = [(end_state, l+1, '')]
    res = ''

    while memory[last][l]:
        _, emit = memory[last][l]
        path.append((last, l, emit))
        state_type = last[0]

        last, emit = memory[last][l]
        res = emit+res
        if state_type in ['M', 'I']:
            l = l-1
        elif state_type == 'SW':
            l = l-2

    path.append((('S', 0), 0, ''))
    path.reverse()

    return dp_score, memory, dp, res, path


def HMM_Viterbi_Learning(seqs, Alphabet, States, Transition, Emission):
    '''
    :param seqs: aligned dataset with same length. Used for HMM training.
    :param Alphabet: dictionary all domain appear in sequence
    :param States: dictionary of all state of HMM
    :param Transition: dictionary of transition rate between two state
    :param Emission: dictionary of emission rate of states
    '''
    # get optimized path of each sequence
    hidden_paths = []
    for seq in seqs:
        dp_score, memory, dp, res, path = align_to_SWAP_HMM(
            seq, Alphabet, States, Transition, Emission)
        hidden_paths.append(path)

    # Adapt tool function from construct_SWAP_HMM function

    def create_state_indices(valid_position):
        '''
        :param valid_position list showing position conservation
        :Given a list showing position is insertion or not, this function return a list of states in HMM
        '''
        index = 0
        state_indices = [('S', index), ('I', index)]
        for p in valid_position:
            if not p:
                continue
            index += 1
            state_indices.append(('M', index))
            state_indices.append(('D', index))
            # (need comment)
            if index >= 2:
                state_indices.append(('SW', index))
            state_indices.append(('I', index))
        index += 1
        state_indices.append(('E', index))
        return state_indices

    def create_state_counts(StateIndices):
        '''
        This function create a frame work illustrating state-state transfer.
        '''
        state_counts = {}
        # Get every pair of states, examine whether there is a one-way transition relationship between them.
        # (s,i) is always regard as the frontal one and (ns,ni) is regard as the posterior one.
        for s, i in StateIndices:
            for ns, ni in StateIndices:
                if s == 'S':
                    if ns == 'I' and ni == 0:
                        state_counts[(s, i), (ns, ni)] = 0
                    elif ns in ['M', 'D'] and ni == 1:
                        state_counts[(s, i), (ns, ni)] = 0
                    # (need comment)
                    elif ns == 'SW' and ni == 2:
                        state_counts[(s, i), (ns, ni)] = 0
                elif s == 'I':
                    # Insertion state is special as it can produce self transition, so the time of HMM chain
                    # wouldn't change, so i should be equal to ni when self transition happen. But for transition
                    # of changing to Deletion and Match state, ni should be equal to i+1.
                    if i == 0:
                        if ns == 'I' and ni == 0:
                            state_counts[(s, i), (ns, ni)] = 0
                        elif ns in ['M', 'D'] and ni == 1:
                            state_counts[(s, i), (ns, ni)] = 0
                        # (need comment)
                        elif ns == 'SW' and ni == 2:
                            state_counts[(s, i), (ns, ni)] = 0
                    else:
                        if (ns == 'I' and ni == i) or (ns in ['M', 'D', 'E'] and ni == i+1):
                            state_counts[(s, i), (ns, ni)] = 0
                        # (need comment)
                        elif ns == 'SW' and ni == i+2:
                            state_counts[(s, i), (ns, ni)] = 0
                elif s in ['M', 'D']:
                    # If current state is one of Matches or Deltion and its going to transfer to state Insertion,
                    # the time of HMM chain wouln't change either.
                    if (ns == 'I' and ni == i) or (ns in ['M', 'D', 'E'] and ni == i+1):
                        state_counts[(s, i), (ns, ni)] = 0
                    # (need comment)*
                    elif ns == 'SW' and ni == i+2:
                        state_counts[(s, i), (ns, ni)] = 0
                elif s == 'SW':
                    # If current state is one of Matches or Deltion and its going to transfer to state Insertion,
                    # the time of HMM chain wouln't change either.
                    if (ns == 'I' and ni == i) or (ns in ['M', 'D', 'E'] and ni == i+1):
                        state_counts[(s, i), (ns, ni)] = 0
                    # (need comment)*
                    elif ns == 'SW' and ni == i+2:
                        state_counts[(s, i), (ns, ni)] = 0
                else:
                    assert(
                        s == 'E'), 'An unexception error occur, please check your input!'
        return state_counts

    def create_emission_counts(StateIndices, Alphabet):
        '''
        This function create a frame work illustrating emission rate of states.
        '''
        emission_counts = {}
        for state in StateIndices:
            for symbol in Alphabet:
                emission_counts[(state, symbol)] = 0
        return emission_counts

    def create_state_frequencies(StateCounts, StateIndices):
        for i in StateCounts.keys():
            StateCounts[i] += m/15  # use global varible m
            # (need comment)* again
            # (need comment)**  explain why the pesudoCount is important for this model and why state SW was specially treated
            # if i[1][0] == 'SW':
            #     StateCounts[i] += m/10
            if i[1][0] == 'I' and i[0][0] == 'I':
                StateCounts[i] -= m/30
        Totals = {i: 0 for i in StateIndices}

        for key, count in StateCounts.items():
            Totals[key[0]] += count

        result = {}
        for key, count in StateCounts.items():
            result[key] = count/Totals[key[0]] if Totals[key[0]] > 0 else 0

        return result

    def create_emission_frequencies(EmissionCounts, StateIndices):
        for i in EmissionCounts.keys():
            if i[0][0] in ['M', 'I']:
                EmissionCounts[i] += k/10000  # use global varible k
        Totals = {i: 0 for i in StateIndices}

        for key, count in EmissionCounts.items():
            Totals[key[0]] += count

        result = {}
        for key, count in EmissionCounts.items():
            result[key] = count/Totals[key[0]] if Totals[key[0]] > 0 else 0

        return result

    StateCounts = create_state_counts(States)
    EmissionCounts = create_emission_counts(States, Alphabet)

    for path in hidden_paths:
        for step in range(1, len(path)):
            prev_state = path[step-1][0]
            current_state = path[step][0]
            emit = path[step][2]
            StateCounts[(prev_state, current_state)] += 1
            if current_state[0] in ['M', 'I']:
                EmissionCounts[current_state, emit] += 1

    Transition = create_state_frequencies(StateCounts, States)
    Emission = create_emission_frequencies(EmissionCounts, States)

    return Transition, Emission


def merge_hidden_paths(pattern_dict):
    '''
    :param paths: a dictionary of pattern to its hidden path which follow the following format: 
        {'BBC':[(('S', 0), 0, ''), (('I', 0), 1, 'B'), (('M', 1), 2, 'B'), (('D', 2), 2, '-'), (('M', 3), 3, 'C'), (('E', 4), 4, '')]}
    :return: a list of aligned sequence ,all of the sequences are of the same length
    '''
    # list of state where insert happen, in format('I', n)

    # paths=[]
    # for k in pattern_dict.keys():
    #     paths.append(pattern_dict(k))
    insert_bin = []
    for pattern, path in pattern_dict.items():
        for state in path:
            if state[0][0] == 'I' and state[0] not in insert_bin:
                insert_bin.append(state[0])

    path_dicts = []
    for pattern, path in pattern_dict.items():
        # create tmp_dict to put normal alignment result of Hidden path and Insertion with their location
        tmp_dict = {'pattern': pattern, 'no_insert_region': ''}
        for state in path:
            if state[0][0] == 'I':
                if state[0] not in tmp_dict:
                    tmp_dict[state[0]] = state[2]
                else:
                    tmp_dict[state[0]] += state[2]
            else:
                tmp_dict['no_insert_region'] += state[2]
        path_dicts.append(tmp_dict)

    # now we want to update alignment result in insert region
    insert_length = {}
    for insert_state in insert_bin:
        # find regions in one insert state in all paths
        insertions = []
        for path_dict in path_dicts:
            if insert_state in path_dict:
                insertions.append(path_dict[insert_state])
            else:
                path_dict[insert_state] = ''

        lengths = [len(x) for x in insertions]
        longest_insertion = insertions[lengths.index(max(lengths))]
        insert_length[insert_state] = lengths
        # print(longest_insertion)

        # Here we make an assumption that the lognest insertion would contain all domain that
        # appear in this insertin state, which is biologically make sense as the region that have
        # vairous domains would usually not be an Instertion state but a Match state.
        for path_dict in path_dicts:
            if path_dict[insert_state] != '':
                insert = path_dict[insert_state]
                new = ''
                j = 0
                for i in range(len(longest_insertion)):
                    if j < len(insert):
                        if insert[j] == longest_insertion[i]:
                            new += longest_insertion[i]
                            j += 1
                        else:
                            new += '-'
                    else:
                        new += '-'*(len(longest_insertion)-j)
                path_dict[insert_state] = new
            else:
                path_dict[insert_state] = '-'*len(longest_insertion)
        # print(path_dicts)

    # now update all sequence
    for path_dict in path_dicts:
        path_dict['aligned'] = path_dict['no_insert_region']

    for insert_state in insert_bin:
        ind = insert_state[1]
        for path_dict in path_dicts:
            old = path_dict['aligned']
            path_dict['aligned'] = old[:ind]+path_dict[insert_state]+old[ind:]

    align_result = {path_dict['pattern']: path_dict['aligned']
                    for path_dict in path_dicts}
    # for r in path_dicts:
    #     print(r)
    # print('-'*50)
    # for r in align_result.items():
    #     print(r)

    return align_result

def set_symbol_for_domains(domain_alignment_result):
    '''
    :param domain_alignment_result: output from greedy MSA alignment
    :return: 
    '''
    patterns = list(map(lambda x: x.split('---'),
                    list(domain_alignment_result.keys())))
    symbol_list = [chr(i)for i in range(65, 91)]
    convert_dict = {}
    convert_dict_rev = {}
    all_domains = []
    for l in patterns:
        all_domains.extend(l)

    all_domains = list(set(all_domains))
    length = len(all_domains)
    for i in range(length):
        convert_dict_rev[symbol_list[i]] = all_domains[i]
        convert_dict[all_domains[i]] = symbol_list[i]

    new_alignment_result = {}
    for i in list(domain_alignment_result.keys()):
        pattern = list((i.split('---')))
        new_pattern = []
        for s in pattern:
            new_pattern.append(convert_dict[s])
        new_alignment_result['---'.join(new_pattern)
                             ] = domain_alignment_result[i]

    return convert_dict_rev, convert_dict, new_alignment_result


def split_to_corresponding_region(align_pattern):
    '''
    split aligned pattern into comparable region
    '''
    values = list(align_pattern.values())
    max_region_number = (len(values[0])*2+1)
    regions_dict = {p: [] for p in list(align_pattern.values())}

    def assign_to_frame(pattern):
        length = len(pattern)
        frame = [''] * (2*length+1)
        for i in range(length):
            frame[2*i+1] = pattern[i]
        left_side = True
        belong_to_last = False
        for i in range(len(frame)):
            if left_side == True:
                if frame[i+1] in ['', '-']:
                    frame[i] = ''
                else:
                    frame[i] = '.'
                    left_side = False
            else:
                if belong_to_last == True:
                    if frame[i] not in ['', '-']:
                        # frame[i] = '.'
                        belong_to_last = False
                    else:
                        frame[i] = '.'
                        belong_to_last = False
                else:
                    if frame[i] not in ['', '-']:
                        # frame[i] = '.'
                        belong_to_last = True
                    else:
                        frame[i] = ''
        return frame

    for pattern in align_pattern.values():
        frame = assign_to_frame(pattern)
        regions_dict[pattern] = frame

    return regions_dict


def parse_domain_info_to_region_info(pattern_dict, regions_dict, domain_alignment_result, align_pattern, all_domains=None):
    '''
    align regions of sequences into their setted position to aligned to other region sequence in that position.
    '''
    def convert_back(id):
        c_id = list(id)
        # for i in id:
        #     c_id.append(convert_dict_rev[i])
        res = '---'.join(c_id)
        return res

    def rm_lower(string):
        newstring = ''.join([i for i in string if not i.islower()])
        return newstring

    def rm_digit(string):
        newstring = ''.join([i for i in string if not i.isdigit()])
        return newstring

    def parse_domain_seq_to_linker(domain_seq, total_length):
        domain_region = list(map(lambda x: x[1:], domain_seq))
        linkers_region = []
        start = 0
        for region in domain_region:
            end = region[0]-1
            linkers_region.append([start, end])
            start = region[1]
        linkers_region.append([start, total_length-1])
        return linkers_region

    original_ids = list(pattern_dict.keys())

    id_switch = {original_id: align_pattern[rm_digit(
        rm_lower(original_id))] for original_id in original_ids}

    sequences_for_HMM = {}
    new_frame = []
    for id in original_ids:
        number = len(domain_alignment_result[convert_back(id)]['seqs'])


        frame = regions_dict[id_switch[id]]
        sequences_for_HMM[id] = {}
        for i in range(number):
            curr_frame_new = []
            total_seq = domain_alignment_result[convert_back(id)]['seqs'][i]
            domain_seq = all_domains[domain_alignment_result[convert_back(id)]['ids'][i]]


            total_length = len(
                domain_alignment_result[convert_back(id)]['seqs'][i])
            linker_region = parse_domain_seq_to_linker(domain_seq, total_length)

            linker_cnt = 0
            domain_cnt = 0
            for f in range(len(frame)):
                if frame[f] == '':
                    curr_frame_new.append('')
                elif frame[f] == '.':
                    lstart, lend = linker_region[linker_cnt]
                    curr_frame_new.append(total_seq[lstart:lend+1])
                    linker_cnt += 1
                elif frame[f] not in ['.', '']:
                    dstart, dend = domain_seq[domain_cnt][1:]
                    curr_frame_new.append(total_seq[dstart:dend+1])
                    domain_cnt += 1
            new_frame.append(curr_frame_new)
    return new_frame



def category_to_HMM_final(domain_alignment_result, all_domains=None, out=None, id_seq_dict=None, mode='real'):
    print('------------------------------------------------')
    print('HMM domain level alignment started...')
    print('------------------------------------------------')


    convert_dict_rev, convert_dict, domain_alignment_result = set_symbol_for_domains(
        domain_alignment_result)

    patterns = list(map(lambda x: x.replace('-', ''),
                    list(domain_alignment_result.keys())))


    def rm_digit(string):
        newstring = ''.join([i for i in string if not i.isdigit()])
        return newstring


    patterns = list(map(rm_digit, patterns))
    Alphabet = sorted(list(set(list(''.join(patterns)))))
    count_of_pattern = [len(domain_alignment_result[key]['ids'])
                        for key in list(domain_alignment_result.keys())]


    most_pattern = patterns[count_of_pattern.index(max(count_of_pattern))]

    States, Transition, Emission = construct_SWAP_HMM(
        0.35, Alphabet, [most_pattern])


    seqs = []
    for i in range(len(patterns)):
        seqs.extend([patterns[i]]*count_of_pattern[i])


    for i in range(100):
        Transition, Emission = HMM_Viterbi_Learning(
            seqs, Alphabet, States, Transition, Emission)

    pattern_dict = {}
    for p in patterns:
        path = align_to_SWAP_HMM(p, Alphabet, States, Transition, Emission)[4]
        pattern_dict[p] = path

    align_pattern ={}
    for p, a in merge_hidden_paths(pattern_dict).items():
        align_pattern[p] = a

    regions_dict = split_to_corresponding_region(align_pattern)

    sequences_for_HMM = parse_domain_info_to_region_info(
        pattern_dict, regions_dict, domain_alignment_result, align_pattern, all_domains)

    res = domain_aware_greedy_MSA(all_domains=None, id_seq_dict=None, mode='hmm', sequences_for_HMM=sequences_for_HMM)
    output = {}
    cnt = 0
    if mode != 'simulate':
        for id in list(all_domains.keys()):
            if id in id_seq_dict:
                output['>|'+id] = ''
                for j in range(len(res)):
                    output['>|'+id] += res[j][cnt][0]
                cnt += 1
    else:
        for id in list(all_domains.keys()):
            output['>|'+id] = ''
            for j in range(len(res)):
                output['>|'+id] += res[j][cnt][0]
            cnt += 1
    print('Outputting result....')
    with open(out+'/result_HMM.fasta', 'w') as outfile:
        for id in output:
            outfile.write(id+'\n')
            outfile.write(output[id]+'\n')
        outfile.close()
    return output
