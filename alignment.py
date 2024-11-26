import os
from Bio import Align
from Bio.Align import substitution_matrices
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

os.chdir(os.path.dirname(os.path.abspath(__file__)))

def read_fasta(file_path):
    ids = []
    sequences = []
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                ids.append(line)
            else:
                sequences.append(line)
    return ids, sequences

#preload data
human_ids, human_sequences = read_fasta("human.fa")
mouse_ids, mouse_sequences = read_fasta("mouse.fa")
example_ids, example_sequences = read_fasta("example.fa")

def calculate_negative_similarity(target_seq, alignments):
    scores = np.zeros(len(target_seq))
    
    for alignment in alignments:
        aligned_target, aligned_subject = alignment[0], alignment[1]
        for i, (t_res, s_res) in enumerate(zip(aligned_target, aligned_subject)):
            if i < len(scores) and t_res != "-" and s_res != "-":
                if t_res == s_res:
                    scores[i] += 1

    max_score = len(alignments)
    negative_similarity = max_score - scores
    
    return negative_similarity

def plot_negative_similarity(target_id, target_seq, negative_similarity):

    positions = list(range(1, len(target_seq) + 1))
    
    plt.figure(figsize=(10, 5))
    plt.bar(positions, negative_similarity, color='skyblue', edgecolor='black')
    plt.xlabel("Sequence position in target")
    plt.ylabel("Negative similarity")
    plt.title("Negative similarity of target sequence to human alignments")
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plot_name = f"{re.sub(r">", "", target_id)}_negative_similarity.png"
    plt.savefig(plot_name)
    plt.show()

def calculate_identity(alignment):
    seg1 = alignment[0]
    seg2 = alignment[1]
    
    total_length = len(seg1) if len(seg1) == len(seg2) else 0
    
    identical_count = sum(1 for i, j in zip(seg1, seg2) if i == j)
        

    return (identical_count / total_length) * 100 if total_length > 0 else 0

def alignment_and_data_creation(target_id, target_seq):
    
    data = []
    
    data.append({"Organism": "Target", "Id": target_id, "Sequence": target_seq})
    
    for h_id, h_seq in zip(human_ids, human_sequences):
        data.append({'Organism': 'Human', 'Id': h_id, 'Sequence': h_seq})

    for m_id, m_seq in zip(mouse_ids, mouse_sequences):
        data.append({'Organism': 'Mouse', 'Id': m_id, 'Sequence': m_seq})

    matrix_choice = input("""
Please choose a substitution matrix:
0. BLOSUM45 (distantly related proteins)
1. BLOSUM62 (midrange)
2. BLOSUM80 (more closely related proteins)\n
""")

    matrix_map = {
        "0": "BLOSUM45",
        "1": "BLOSUM62",
        "2": "BLOSUM80"
    }
    
    selected_matrix = matrix_map.get(matrix_choice, "BLOSUM62")
    
    
    alignment_choice = input("""
Please choose an alignment mode:
0. Automatic (biopython chooses appropriate mode)
1. Local alignment
2. Global alignment\n
""")
    
    alignment_map = {
        "0": None,
        "1": "local",
        "2": "global"
    }
    
    selected_alignment = alignment_map.get(alignment_choice, None)
    
    aligner = Align.PairwiseAligner()
    
    if selected_alignment:
        aligner.mode = selected_alignment
        
    aligner.substitution_matrix = substitution_matrices.load(selected_matrix)

    for entry in data[1:]:
        alignment = aligner.align(target_seq, entry['Sequence'])[0]
        
        entry['Pairwise Alignment'] = [alignment[0], alignment[1]]
        entry['Alignment_score'] = alignment.score
        entry['Identity_percentage'] = calculate_identity(alignment)

    df = pd.DataFrame(data)

    print(df)
    print(f"Substitution matrix: {selected_matrix}\nAlingment mode: {aligner.mode}")
    
    while True:
        try:
            user_action = int(input("""
Please choose an action:
1. Print summary statistics
2. Save DataFrame to CSV
3. Graph negative similarity with human database
4. Exit\n
"""))
            
            if user_action == 1:
                calculate_summary_statistics(df)
            elif user_action == 2:
                csv_file_path = "alignment_results.csv"
                df.to_csv(csv_file_path, index=False)
                print(f"Data saved to {csv_file_path}")
            elif user_action == 3:
                human_alignments = df[df["Organism"] == "Human"]["Pairwise Alignment"]
                negative_similarity = calculate_negative_similarity(target_seq, human_alignments)
                plot_negative_similarity(target_id, target_seq, negative_similarity)
            elif user_action == 4:
                print("Exiting...")
                break
            else:
                print("Invalid option. Please select a number between 1 and 3.")
        except ValueError:
            print("Input is not a valid integer. Please try again.")
    
    return df

def custom_sequence_handle():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    print(f"This is your current directory: {current_dir}")
    
    file_path = input("Please input the path to your fasta file: ")
    
    if os.path.isfile(file_path):
        print(f"File '{file_path}' found and ready to process.")
        ids, sequences = read_fasta(file_path)
        target_ids = ids[0]
        target_seq = sequences[0]
        return alignment_and_data_creation(target_ids, target_seq)
    else:
        print(f"File '{file_path}' not found. Please try again.")
        return custom_sequence_handle()

def handle_examples(option):
    print(f"You selected option: {option}, with sequence:\n{example_sequences[option]}")
    target_id = example_ids[option]
    target_seq = example_sequences[option]
    return alignment_and_data_creation(target_id, target_seq)

def calculate_summary_statistics(df):
    
    summary_stats = {}
    
    grouped = df.groupby('Organism', group_keys=True)['Identity_percentage']
    
    for organism, group in grouped:
        if not organism == "Target":

            max_identity = group.max()
            min_identity = group.min()
            mean_identity = group.mean()
            median_identity = group.median()
            
            summary_stats[organism] = {
                'max_identity': max_identity,
                'min_identity': min_identity,
                'mean_identity': mean_identity,
                'median_identity': median_identity
            }
        
            print(f"\nStatistics for {organism}:")
            print(f"  Max Identity: {max_identity:.2f}%")
            print(f"  Min Identity: {min_identity:.2f}%")
            print(f"  Mean Identity: {mean_identity:.2f}%")
            print(f"  Median Identity: {median_identity:.2f}%")

    max_organism = max(summary_stats.keys(), key = lambda org: abs(summary_stats[org]['max_identity']))
    mean_organism = max(summary_stats.keys(), key = lambda org: abs(summary_stats[org]["mean_identity"]))
        
    if max_organism == mean_organism:
        print(f"\nThe target sequence most resembles {max_organism.lower()}")
    else:
        print(f"The maximum sequence identity indicates that the target sequence most resembles {max_organism.lower()}")
        print(f"The mean sequence identity indicates that the target sequence most resembles {mean_organism.lower()}")



def menu():
    while True:
        try:
            ans = int(input(
                """
Please choose one of the following options:
0. Custom sequence
1. Ocrelizumab Variable heavy chain
2. Infliximab variable heavy chain
3. Muromonab variable heavy chain
4. Bevacizumab variable heavy chain
5. Caplacizumab variable heavy chain\n
"""
            ))
            if ans == 0:
                return custom_sequence_handle()
            elif 1 <= ans <= 5:
                return handle_examples(ans - 1)
            else:
                print("Invalid option. Please select a number between 0 and 5.")
        except ValueError:
            print("Input is not a valid integer. Please try again.")


if __name__ == "__main__":
    menu()
