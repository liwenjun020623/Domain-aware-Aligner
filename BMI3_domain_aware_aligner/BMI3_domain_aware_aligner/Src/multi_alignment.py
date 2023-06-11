"""
BMI3 Group5 Project4: Domain-aware aligner
"""


import numpy as np
from math import inf
import blosum as bl
from copy import deepcopy
from Src.uniprot_domain_processor import *


"""
PARSE INPUT
"""


def parse_txt_to_dict(path):
    """
    :param path: str path of a txt file of fasta format
    :return: dict {id: seq, ...}
    """
    with open(path, 'r') as infile:
        lines = infile.readlines()
    fasta_dict = {}
    curr_seq = ''
    curr_id = None
    for row in lines:
        if row[0] == '>':
            # record, except for beginning
            if curr_id is not None:
                fasta_dict[curr_id] = curr_seq
            # next seq
            curr_seq = ''
            curr_id = row.split('|')[1]
        else:
            curr_seq += row.strip().replace('-', '')
    # add last
    fasta_dict[curr_id] = curr_seq
    return fasta_dict


"""
CORE ALGORITHM: PAIR-WISE GLOBAL ALIGNMENT WITH AFFINE GAP PENALTY
"""


def affine_gap_penalties_pair_wise_alignment(seq1, seq2=None, sigma=11, epsilon=1,
                                             mode='pairwise', alignment=None, outmode='score+align'):
    """
    :param mode: str select between 'pairwise' & 'profile' modes
    :param seq1: str input protein sequence 1
    :param seq2: str input protein sequence 2
    :param sigma: int penalty for opening a gap
    :param epsilon: int penalty for extending a gap
    :param alignment: list of alignment strings
    :return: tuple (align1: str alignment for seq1, align2: str alignment for seq2, score: float best alignment score)
    """
    def gen_backtrack(seq1, seq2, sigma=11, epsilon=1, mode='pairwise', profile=None):
        """
        :param profile: 2d nested list of profile matrix of previous alignment
        :param mode: str select between 'pairwise' & 'profile' modes
        :param seq1: str input protein sequence 1
        :param seq2: str input protein sequence 2
        :param sigma: int penalty for opening a gap
        :param epsilon: int penalty for extending a gap
        :return: tuple (backtracks: list [3 np.array backtracks for gap1, match, and gap2], score: float best alignment score)
        """
        def weighted_score_with_profile(profile_column1, profile_column2, penalty):
            """
            :param profile_column1: list of float the column of alignment profile matrix 1 at corresponding position
            :param profile_column2: list of float the column of alignment profile matrix 2 at corresponding position
            :param penalty: dictionary penalty matrix
            :return: weighted_score: float score of weighted score based on profile
            """
            amino_arr = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                         'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
            weighted_score = 0
            if len(profile_column1) == 1:
                for row in range(len(profile_column1)):
                    weighted_score += penalty[profile_column1 + amino_arr[row]] * profile_column2[row]
            else:
                for row1 in range(len(profile_column1)):
                    for row2 in range(len(profile_column2)):
                        weighted_score += penalty[amino_arr[row1] + amino_arr[row2]] * profile_column2[row2] * profile_column1[row1]
            return weighted_score
        # init penalty
        penalty = bl.BLOSUM(62)
        len1 = len(seq1)
        if mode == 'pairwise':
            len2 = len(seq2)
        elif mode == 'profile':
            assert profile is not None, 'profile parameter required under profile mode!'
            len2 = len(profile[0])
        elif mode == 'profile2':
            assert profile is not None, 'profile parameter required under profile2 mode!'
            len1 = len(profile[0][0])
            len2 = len(profile[1][0])
        else:
            assert False, 'Invalid mode! Please select from "pairwise" and "profile"!'
        # init scores
        match_scores = -inf * np.ones((len1+1, len2+1))
        gap1_scores = -inf * np.ones((len1+1, len2+1))
        gap2_scores = -inf * np.ones((len1+1, len2+1))
        # init backtracks
        # what are backtracks:
        # they are matrices of the same dimension as score matrices, used to record the best previous position
        match_back = [[0 for _ in range(len2+1)] for _ in range(len1+1)]
        # first row of match_back should be 2, representing gap2
        match_back[0] = [2 for _ in range(len2 + 1)]
        gap1_back = [[0 for _ in range(len2+1)] for _ in range(len1+1)]
        gap2_back = [[0 for _ in range(len2+1)] for _ in range(len1+1)]
        # at start position, all scores are 0
        gap1_scores[0, 0] = 0
        gap2_scores[0, 0] = 0
        match_scores[0, 0] = 0
        # init first gap openings
        gap1_scores[1, 0] = -sigma
        gap2_scores[0, 1] = -sigma
        match_scores[1, 0] = -sigma
        match_scores[0, 1] = -sigma
        # init lateral gap extensions: technically, we have 3 score matrices, we first calculate along all edges
        for pos1 in range(2, len1 + 1):
            gap1_scores[pos1, 0] = gap1_scores[pos1 - 1, 0] - epsilon
            match_scores[pos1, 0] = gap1_scores[pos1, 0]
        for pos2 in range(2, len2 + 1):
            gap2_scores[0, pos2] = gap2_scores[0, pos2 - 1] - epsilon
            match_scores[0, pos2] = gap2_scores[0, pos2]
        # start dynamic traversing: technically, we start to dynamically fill in the 3 matrices
        # unlike normal alignment, for affine gap alignment, we need 3 matrices,
        # so for each point, we consider 3 states: in gap1, match, in gap2
        for pos1 in range(1, len1 + 1):
            for pos2 in range(1, len2 + 1):
                # First, consider we are now in gap1
                # to best get to gap 1, do we open new gap? or extend existing gap?
                # if we extend, we move horizontally in gap1 state
                # if we open, we move diagonally from match state to gap1 state
                gap1_extend_score = gap1_scores[pos1 - 1, pos2] - epsilon
                gap1_open_score = match_scores[pos1 - 1, pos2] - sigma
                gap1_best_score = np.max([gap1_open_score, gap1_extend_score])
                gap1_scores[pos1, pos2] = gap1_best_score
                # if we open gap:
                if gap1_best_score == gap1_open_score:
                    # we record at gap1 backtrack as from match
                    gap1_back[pos1][pos2] = 1
                else:
                    # if extend
                    gap1_back[pos1][pos2] = 0
                # Then, consider we are now in gap2
                # to best get to gap 2, do we open new gap? or extend existing gap?
                # if we extend, we move horizontally in gap2 state
                # if we open, we move diagonally from match state to gap2 state
                gap2_extend_score = gap2_scores[pos1, pos2 - 1] - epsilon
                gap2_open_score = match_scores[pos1, pos2 - 1] - sigma
                gap2_best_score = np.max([gap2_open_score, gap2_extend_score])
                gap2_scores[pos1, pos2] = gap2_best_score
                # if we should open gap:
                if gap2_best_score == gap2_open_score:
                    # we record at gap1 backtrack as from match
                    gap2_back[pos1][pos2] = 1
                else:
                    gap2_back[pos1][pos2] = 2
                # Last, consider we are now in match state
                # do we continue with gaps / do a match?
                match_match_score = None
                if mode == 'pairwise':
                    match_match_score = match_scores[pos1 - 1, pos2 - 1] + penalty[seq1[pos1-1] + seq2[pos2-1]]
                elif mode == 'profile':
                    # get the column to calculate weighted score from profile
                    profile_col = [row[pos2-1] for row in profile]
                    match_match_score = match_scores[pos1 - 1, pos2 - 1] + weighted_score_with_profile(seq1[pos1-1], profile_col, penalty)
                elif mode == 'profile2':
                    profile_col1 = [row[pos1-1] for row in profile1]
                    profile_col2 = [row[pos2-1] for row in profile2]
                    match_match_score = match_scores[pos1 - 1, pos2 - 1] + weighted_score_with_profile(profile_col1, profile_col2, penalty)
                best_score = np.max([gap1_best_score,
                                     match_match_score,
                                     gap2_best_score])
                match_scores[pos1, pos2] = best_score
                if best_score == match_match_score:
                    match_back[pos1][pos2] = 1
                elif best_score == gap2_best_score:
                    match_back[pos1][pos2] = 2
                elif best_score == gap1_best_score:
                    match_back[pos1][pos2] = 0
        return [gap1_back, match_back, gap2_back], match_scores[len1, len2]

    def gen_alignment_from_backtrack(backtracks, seq1, seq2, mode='pairwise', profile=None):
        """
        :param profile: 2d nested list of profile matrix of previous alignment
        :param mode: str select between 'pairwise' & 'profile' modes
        :param backtracks: list [3 np.array backtracks for gap1, match, and gap2]
        :param seq1: str input protein sequence 1
        :param seq2: str input protein sequence 2
        :return: tuple (align1: str alignment for seq1, align2: str alignment for seq2)
        """
        # start at the end of backtrack
        pos1 = len(seq1)
        if mode == 'pairwise':
            pos2 = len(seq2)
        elif mode == 'profile':
            assert profile is not None, 'profile parameter required under profile mode!'
            pos2 = len(profile[0])
            # if we align to profile, seq2 will be a dummy sequence here just to record the gaps
            seq2 = ''.join(['X' for _ in range(len(profile[0]))])
        elif mode == 'profile2':
            assert profile is not None, 'profile parameter required under profile2 mode!'
            profile1, profile2 = profile
            pos1 = len(profile1[0])
            pos2 = len(profile2[0])
            seq1 = ''.join(['X' for _ in range(len(profile1[0]))])
            seq2 = ''.join(['X' for _ in range(len(profile2[0]))])
        else:
            assert False, 'Invalid mode! Please select from "pairwise" and "profile"!'
        # extract backtracts
        gap1_back, match_back, gap2_back = backtracks
        # init alignment strings
        align1 = ''
        align2 = ''
        state = 1  # match
        # while not finished, we go back one step at a time
        while pos1 > 0 or pos2 > 0:
            # if we now in gap1
            if state == 0:
                # if this is a gap opening, we go back to match state
                if gap1_back[pos1][pos2] == 1:
                    state = 1
                # we move back one step in seq1, but not seq2 because we are in gap1
                pos1 -= 1
                align1 += seq1[pos1]
                align2 += '-'
            # if we now in gap2
            elif state == 2:
                # if this is a gap2 opening, we go back to match state
                if gap2_back[pos1][pos2] == 1:
                    state = 1
                # we move back one step in seq2, but not seq1 because we are in gap2
                pos2 -= 1
                align1 += '-'
                align2 += seq2[pos2]
            # if we now in match state
            elif state == 1:
                # what did we do last time?
                # did we come from match / gap1 closing / gap2 closing?
                prev_state = match_back[pos1][pos2]
                # if we came from a match,
                # we go back one step in BOTH seq1 and seq2
                if prev_state == 1:
                    pos1 -= 1
                    pos2 -= 1
                    align1 += seq1[pos1]
                    align2 += seq2[pos2]
                # if we came from either gap,
                # we go one step back to gap1 / gap2 state
                elif prev_state in [0, 2]:
                    state = prev_state
        # when we are at the start, we return results
        return align1, align2

    def alignment_to_profile(alignment):
        """
        :param alignment: list of alignment strings
        :return: matrix: 2d-nested list of profile matrix
        """
        matrix = []
        for base in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                     'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']:
            mat = []
            for pos in range(len(alignment[0])):
                col = [row[pos] for row in alignment]
                mat.append(col.count(base) / len(alignment))
            matrix.append(mat)
        return matrix
    # initialize results
    align1 = align2 = score = None
    assert mode in ['pairwise', 'profile', 'profile2'], \
        'Invalid Mode: available modes: [pairwise, profile, profile2]'
    profile = None
    if mode == 'pairwise':
        # normal pairwise procedure
        profile = None
    if mode == 'profile':
        assert alignment is not None, 'alignment parameter required under "profile" mode!'
        # first generate profile of the alignment
        profile = alignment_to_profile(alignment)
    if mode == 'profile2':
        assert alignment is not None, 'alignment parameter required under "profile2" mode!'
        align1, align2 = alignment
        profile1 = alignment_to_profile(align1)
        profile2 = alignment_to_profile(align2)
        profile = [profile1, profile2]
    backtracks, score = gen_backtrack(seq1, seq2, sigma, epsilon, mode, profile)
    # initialize align1 & 2 as empty string and only really compute then under score+align mode to speed up
    align1 = align2 = ''
    if outmode == 'score+align':
        align1, align2 = gen_alignment_from_backtrack(backtracks, seq1, seq2, mode, profile)
    return align1[::-1], align2[::-1], score


"""
MAIN ALGORITHM: GREEDY MSA DOMAIN AWARE ALIGNER (with built-in logger)
"""


# Tree Data Structure used in Hierarchical Clustering
class Tree:
    def __init__(self):
        self.root = None
        self.forward_dict = {}
        self.reverse_dict = {}
        self.adj_dict = {}
        self.nodes = set()
        self.node_weight_dict = {}

    def addConnection(self, node1, node2, dis):
        if node1 not in self.adj_dict:
            self.adj_dict[node1] = {node2: dis}
        else:
            self.adj_dict[node1][node2] = dis
        if node2 not in self.adj_dict:
            self.adj_dict[node2] = {node1: dis}
        else:
            self.adj_dict[node2][node1] = dis
        self.nodes.add(node1)
        self.nodes.add(node2)

        if node1 not in self.forward_dict:
            self.forward_dict[node1] = {node2: dis}
        else:
            self.forward_dict[node1][node2] = dis
        if node2 not in self.reverse_dict:
            self.reverse_dict[node2] = {node1: dis}
        else:
            self.reverse_dict[node2][node1] = dis

    def getConnections(self, node):
        return self.adj_dict[node]

    def getFather(self, node):
        return self.reverse_dict[node]

    def getChildren(self, node):
        return self.forward_dict[node]

    def addNode(self, node, weight=0):
        self.node_weight_dict[node] = weight
        self.adj_dict[node] = {}
        self.nodes.add(node)

    def getNodeCount(self):
        return len(self.nodes)


def domain_aware_greedy_MSA(all_domains, id_seq_dict, sigma=11, epsilon=1, mode='seq', sequences_for_HMM=None):
    """
    :param all_domains: dict pre-constructed dictionary from UniProt
    :param id_seq_dict: dict {id: seq, ...}
    :param sigma: int penalty for opening a gap
    :param epsilon: int penalty for extending a gap
    :return: dict {structure_identifier: {}, ...}
    """
    def categorize_seqs_by_domain_info(all_domains, id_seq_dict):
        """
        :param all_domains: dict preconstructed dictionary from UniProt
        :param id_seq_dict: dict {UniProt ID: corresponding sequence}
        :return: dict {strucure_identifier: {'seqs': [sequences of this structure], 'ids': [ids of sequences]}, ...}
        """
        # categorize seqs by domain structures
        categories = {}
        # extract ids & seqs from input dictionary
        ids = list(id_seq_dict.keys())
        seqs = list(id_seq_dict.values())
        for i in range(len(ids)):
            curr_id = ids[i]
            curr_seq = seqs[i]
            _, domain_structure_list = sequence_to_domain_structure(curr_id, curr_seq, all_domains)
            # domain_structure_list = curr_domain_info['structure_list_']
            domain_structure_list_no_linkers = [name for name in domain_structure_list if name[-1] != '_']
            # unique identifier of a sequence is all domain names separated by '---' Example: 'KH 1---DH---VHS' 3 domains
            if not domain_structure_list_no_linkers:
                structure_identifier = '_unknown_'
            else:
                structure_identifier = '---'.join(domain_structure_list_no_linkers)
            if structure_identifier not in categories:
                categories[structure_identifier] = {'seqs': [curr_seq], 'ids': [curr_id]}
            else:
                categories[structure_identifier]['seqs'].append(curr_seq)
                categories[structure_identifier]['ids'].append(curr_id)
        return categories

    def split_seq_by_structure_dict(structure_arr, seq):
        """
        :param structure_arr: list generated by sequence_to_domain_structure function
        :param seq: sequence to be split
        :return: list of split domains & linkers (linker could be empty string)
        """
        # init res
        split_seq_list = []
        # iterate through each domain & linker
        for structure in structure_arr:
            # start & end are 0-index, remember to use end+1 for indexing
            start, end = structure[1:]
            split_seq_list.append(seq[start: end + 1])
        return split_seq_list

    def update_alignment(alignment, align2_dummy, mode='string'):
        """
        :param mode: str "domain": input list, or "string": input string
        :param alignment: list of old alignment strings
        :param align2_dummy: str align2 generated recording information about new gaps to be added in old alignments
        :return: alignment_new: list of updated alignment strings with new gaps added
        """
        assert mode in ['domain', 'string'], 'Invalid Mode! mode in update alignment must be domain or string!'
        alignment_new = None
        # base case: if input is a single string, not split by domains
        if mode == 'string':
            # init result as 1d list
            alignment_new = []
            for align in alignment:
                align_new = ''
                align_pos = 0
                for pos in range(len(align2_dummy)):
                    if align2_dummy[pos] == '-':
                        align_new += '-'
                    else:
                        align_new += align[align_pos]
                        align_pos += 1
                alignment_new.append(align_new)
        # if input is list split by domains, we do string mode for each fragment
        elif mode == 'domain':
            # init result as 2d list
            alignment_new = [[] for _ in range(len(alignment))]
            for frag_pos in range(len(alignment[0])):
                # get a column in alignment
                curr_align_ = [align[frag_pos] for align in alignment]
                curr_dummy = align2_dummy[frag_pos]
                # for each fragment, we can use string mode to get a column
                new_align = update_alignment(curr_align_, curr_dummy, mode='string')
                # now we have 1 column, but we need to append by row
                for i in range(len(new_align)):
                    alignment_new[i].append(new_align[i])
        return alignment_new

    def hierarchical_clustering(D_dict, n):
        """
        :param D_dict: distance dictionary for all pairwise alignments
        :param n: number of sequences
        :return: nested list trace, recording clusters to be merged in each step to guide alignment merging
        """
        trace = []

        def average_cluster_dis(D_dict, cluster1, cluster2):
            """
            :param D_dict: distance dictionary for all pairwise alignments
            :param cluster1: list cluster1 nodes
            :param cluster2: list cluster2 nodes
            :return: float distance between 2 clusters
            """
            dis = 0
            for node1 in cluster1:
                for node2 in cluster2:
                    dis += D_dict[node1][node2]
            return dis / (len(cluster1) * len(cluster2))

        def closest_clusters(D_copy_full, clusters):
            """
            :param D_copy_full: distance dictionary for all pairwise alignments
            :param clusters: nested list all clusters
            :return: two closest clusters
            """
            dis_dict = {key: {key: 0 for key in clusters} for key in clusters}
            # compute all distances
            for key1 in dis_dict:
                for key2 in dis_dict[key1]:
                    dis_dict[key1][key2] = average_cluster_dis(D_copy_full, clusters[key1], clusters[key2])
            # find min distance
            best_dis = inf
            best_cluster1 = None
            best_cluster2 = None
            for key1 in dis_dict:
                for key2 in dis_dict[key1]:
                    if dis_dict[key1][key2] < best_dis and key1 != key2:
                        best_dis = dis_dict[key1][key2]
                        best_cluster1 = key1
                        best_cluster2 = key2
            return best_cluster1, best_cluster2, best_dis

        # initialize tree
        tree = Tree()
        for i in range(1, n + 1):
            tree.addNode(i)
        # initialize clusters
        clusters = {i: [i] for i in range(1, n + 1)}
        D_copy_full = deepcopy(D_dict)
        while len(clusters) > 1:
            closest_cluster1, closest_cluster2, closest_dis = closest_clusters(D_copy_full, clusters)
            cluster_new = clusters[closest_cluster1] + clusters[closest_cluster2]
            trace.append([clusters[closest_cluster1], clusters[closest_cluster2]])
            new_node_label = tree.getNodeCount() + 1
            tree.addNode(new_node_label)
            tree.addConnection(new_node_label, closest_cluster1, 0.5 * closest_dis)
            tree.addConnection(new_node_label, closest_cluster2, 0.5 * closest_dis)
            # add Cnew to Clusters
            clusters[new_node_label] = cluster_new
            # add a row/column to D for Cnew by computing D(Cnew, C) for each C in Clusters
            D_dict[new_node_label] = {key: 0 for key in D_dict}
            for key in D_dict:
                current_dis = average_cluster_dis(D_copy_full, clusters[new_node_label], clusters[key])
                D_dict[new_node_label][key] = current_dis
                D_dict[key][new_node_label] = current_dis
            # remove Ci and Cj from Clusters
            clusters.pop(closest_cluster1, None)
            clusters.pop(closest_cluster2, None)
            # remove two rows & cols from D
            for i in [closest_cluster1, closest_cluster2]:
                D_dict.pop(i, None)
                for key in D_dict:
                    D_dict[key].pop(i, None)
        return tree, trace

    def calculate_pairwise_distance_matrix(curr_split_seqs, sigma, epsilon, mode='seq'):
        """
        :param curr_split_seqs: nested list of sequences split by domain composition
        :param sigma: int penalty for opening a gap
        :param epsilon: int penalty for extending a gap
        :return: hamming distance matrix of all pairwise alignments
        """
        def hammingDistance(str1, str2):
            """
            :param str1: seq1
            :param str2: seq1
            :return: int hamming distance between sequences
            """
            dis = 0
            for i in range(0, len(str1)):
                if str1[i] != str2[i]:
                    dis += 1
            return dis
        D_dict = {}
        # for all sequence pairs
        for i in range(len(curr_split_seqs)):
            if mode != 'hmm':
                print('Calculating All Pairwise Distances to ' + str(i+1) + '/'+str(len(curr_split_seqs))+' sequence...')
            for j in range(i+1, len(curr_split_seqs)):
                # first calculate alignment
                curr_seq1 = curr_split_seqs[i]
                curr_seq2 = curr_split_seqs[j]
                align1 = ''
                align2 = ''
                for frag_i in range(len(curr_seq1)):
                    frag_align1, frag_align2, _ = affine_gap_penalties_pair_wise_alignment(curr_seq1[frag_i], curr_seq2[frag_i],
                                                                                           sigma, epsilon, mode='pairwise',
                                                                                           outmode='score+align')
                    align1 += frag_align1
                    align2 += frag_align2
                # then calculate distance
                distance = hammingDistance(align1, align2)
                # record distance in 2 directions
                if i+1 not in D_dict:
                    D_dict[i+1] = {j+1: distance}
                else:
                    D_dict[i+1][j+1] = distance
                if j+1 not in D_dict:
                    D_dict[j+1] = {i+1: distance}
                else:
                    D_dict[j+1][i+1] = distance
        # define distance with onself as 0
        for key in D_dict:
            D_dict[key][key] = 0
        return D_dict

    def merge_two_sequence_clusters(trace, curr_align, mode='seq'):
        """
        :param trace: nested list of clusters to be merged at each step
        :param curr_align: current alignment
        :return: new alignment after merging two sequence clusters
        """
        # get two clusters to be merged
        cluster1, cluster2 = trace
        # get split sequences in those two clusters
        align1 = [curr_align[i - 1] for i in cluster1]
        align2 = [curr_align[j - 1] for j in cluster2]
        # init dummy sequences to record gap insertions
        final_alignment_align_dummy1 = []
        final_alignment_align_dummy2 = []
        # perform pairwise alignment (profile-profile) for on two alignments for each domain / linker
        for frag_i in range(len(align1[0])):
            align1_frag = [align[frag_i] for align in align1]
            align2_frag = [align[frag_i] for align in align2]
            alignment_frag_align_dummy1, alignment_frag_align_dummy2, _ = affine_gap_penalties_pair_wise_alignment(
                ''.join('X' for _ in range(len(align1_frag[0]))),
                ''.join('X' for _ in range(len(align2_frag[0]))),
                sigma, epsilon,
                mode='profile2', alignment=[align1_frag, align2_frag], outmode='score+align')
            final_alignment_align_dummy1.append(alignment_frag_align_dummy1)
            final_alignment_align_dummy2.append(alignment_frag_align_dummy2)
        # update alignment by inserting gaps
        alignment1_new = update_alignment(align1, final_alignment_align_dummy1, mode='domain')
        alignment2_new = update_alignment(align2, final_alignment_align_dummy2, mode='domain')
        ######################
        if mode != 'hmm':
            print('Merged Cluster: ' + str(cluster1) + ' and ' + str(cluster2))
            print('-------------------------------------------------------')
        ######################
        for i in range(len(cluster1)):
            curr_align[cluster1[i]-1] = alignment1_new[i]
        for j in range(len(cluster2)):
            curr_align[cluster2[j]-1] = alignment2_new[j]
        return curr_align
    if mode == 'seq':
        ######################
        print('-------------------------------------------------------')
        print('Domain-based Greedy MSA started...')
        print('-------------------------------------------------------')
        print('Categorizing '+str(len(id_seq_dict))+' sequences')
        ######################
        # split seqs into categories
        categories = categorize_seqs_by_domain_info(all_domains, id_seq_dict)
        ######################
        print('Sequence Categorization Complete')
        print('-------------------------------------------------------')
        print(str(len(categories))+' distinct domain-combinations found, they are (separated by "---"): ')
        for key in list(categories.keys()):
            print(key+' : '+str(len(categories[key]['seqs']))+' sequences')
        print('-------------------------------------------------------')
        ######################
        # init alignment
        alignments = categories.copy()
        # greedy MSA for each categories
        for structure_identifier in categories:
            ######################
            print('Performing greedy MSA on '+str(len(categories[structure_identifier]['seqs']))+' sequences of structure '
                  + structure_identifier)
            print('Calculating Pairwise Distance Matrix, this may take a while...')
            print('-------------------------------------------------------')
            ######################
            # extract seqs & ids for this category
            curr_seqs = categories[structure_identifier]['seqs']
            curr_ids = categories[structure_identifier]['ids']
            # split all seqs by domains
            curr_split_seqs = []
            for i in range(len(curr_seqs)):
                curr_split, _ = sequence_to_domain_structure(curr_ids[i], curr_seqs[i], all_domains)
                curr_split_seqs.append(split_seq_by_structure_dict(curr_split, curr_seqs[i]))
            # calculate distance matrix
            D_dict = calculate_pairwise_distance_matrix(curr_split_seqs, sigma, epsilon)
            # perform hierarchical clustering
            tree, traces = hierarchical_clustering(D_dict, len(curr_split_seqs))
            curr_align = curr_split_seqs
            for trace in traces:
                curr_align = merge_two_sequence_clusters(trace, curr_align)
            # now the alignment process finished, our curr_align is still fragmented (2d list instead of 1d)
            # IF WE WANT DIRECT OUTPUT, we concatenate fragmented alignment right now, but this cost domain information
            concat_final_alignment = [''.join(align_list) for align_list in curr_align]
            # record alignment for this category
            alignments[structure_identifier]['category_alignment'] = concat_final_alignment
            alignments[structure_identifier]['unconcat_align'] = curr_align
    elif mode == 'hmm':
        final_res = []
        curr_seqs = sequences_for_HMM
        gap_info = []
        # for key1 in sequences_for_HMM:
        #     for key2 in sequences_for_HMM[key1]:
        #         curr_seqs.append(sequences_for_HMM[key1][key2])
        for frag_i in range(len(curr_seqs[0])):
            curr_gap = []
            # curr_col = [[curr_seqs[row][frag_i]] for row in range(len(curr_seqs)) if curr_seqs[row][frag_i] != '']
            curr_col_raw = [[curr_seqs[row][frag_i]] for row in range(len(curr_seqs))]
            curr_col = []
            for seq in curr_col_raw:
                if seq == ['']:
                    curr_gap.append(True)
                else:
                    curr_gap.append(False)
                    curr_col.append(seq)
            # curr_col = [seq for seq in curr_col_raw if seq != ['']]
            gap_info.append(curr_gap)
            # curr_seqs = sequences_for_HMM
            curr_split_seqs = curr_col
            # calculate distance matrix
            D_dict = calculate_pairwise_distance_matrix(curr_split_seqs, sigma, epsilon, mode=mode)
            # perform hierarchical clustering
            tree, traces = hierarchical_clustering(D_dict, len(curr_split_seqs))
            curr_align = curr_split_seqs
            for trace in traces:
                curr_align = merge_two_sequence_clusters(trace, curr_align, mode=mode)
            new_curr_align = []
            seq_cnt = 0
            for i in range(len(curr_gap)):
                if curr_gap[i]:
                    new_curr_align.append(['-'*len(curr_align[0][0])])
                else:
                    new_curr_align.append(curr_align[seq_cnt])
                    seq_cnt += 1

            final_res.append(new_curr_align)
        alignments = final_res
    ######################
    if mode != 'hmm':
        print('-------------------------------------------------------')
        print('Domain-based Greedy MSA finished!')
        print('-------------------------------------------------------')
    ######################
    return alignments

