import argparse
import numpy as np
import itertools
import blosum as bl
matrix = bl.BLOSUM(62)

parser = argparse.ArgumentParser(description='This is a sum-of-pairs scoring system.')
parser.add_argument('-i', required=True, type=str, help="the path to read in result fasta file")
args = parser.parse_args()

def readin_resultFile(i:str):
  with open(i) as fl:
    result_seq = []
    content = {}
    entry_id = fl.readline().strip().split('|')[1]
    seq = ''
    for i in fl:
      if '>' in i:
        content[entry_id] = seq
        entry_id = i.strip().split('|')[1]
        seq = ''
      else:
        seq += i.strip().replace('-','*')
    content[entry_id] = seq
    for j in content.values():
      result_seq.append(list(j))
    bench_data = np.array(result_seq)
  return bench_data

def Sum_of_Pairs(bench_data):
  col_num = len(bench_data[0])
  SP_value = 0
  for i in range(1,col_num+1):
    current_col = bench_data[:,i-1:i]
    current_col_list = current_col.tolist()
    current_col_str = []
    for j in current_col_list:
      current_col_str.append(j[0])
    current_combinations = list(itertools.combinations(current_col_str,2))
    current_combinations_str = []
    for item in current_combinations:
      combination_item = item[0] + item[1]
      current_combinations_str.append(combination_item)
    current_col_value = 0
    for combination_item in current_combinations_str:
      current_col_value += matrix[combination_item]
    SP_value += current_col_value
  return SP_value

bench_data = readin_resultFile(args.i)
print("The sum-of-pairs score of this alignment is " + str(Sum_of_Pairs(bench_data)))