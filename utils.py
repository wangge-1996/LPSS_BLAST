from Bio import SeqIO
import json
import pandas as pd
import os

# Extract sequence based on id
def extract_sequences(fasta_file, ids_to_extract):
    extracted_sequences = []
    with open(fasta_file, 'r') as fasta:
        for record in SeqIO.parse(fasta, 'fasta'):
            if record.id in ids_to_extract:
                extracted_sequences.append(record)
    return extracted_sequences

# Save the sequence as a json file format, dedicated to alphafold3
def save_json(sequence, path):
    data = {
        "name": str(sequence.id),
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence":str(sequence.seq)
                }
            }
        ],
        "modelSeeds": [1],
        "dialect": "alphafold3",
        "version": 1
    }
    
    with open(path, 'w') as json_file:
        json.dump(data, json_file, indent=4)

# Get sequence length
def get_length_file(fasta_file):
    extracted_sequences=[]
    with open(fasta_file, 'r', encoding='utf-8') as fasta:
        for record in SeqIO.parse(fasta, 'fasta'):
            extracted_sequences.append([record.id,len(record.seq)])
    
    df = pd.DataFrame(extracted_sequences)
    df.columns = ['reference_id', 'reference_length']
    
    return df

# 定义一个函数来合并区间并计算覆盖范围的长度
def merge_intervals(intervals):
    # 将区间按照起始位置排序
    intervals = sorted(intervals, key=lambda x: x[0])
    # 合并重叠的区间
    merged = []
    for start, end in intervals:
        if not merged or merged[-1][1] < start:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    # 计算合并后区间的总长度
    return sum(end - start for start, end in merged)

# Process blast files
def process_blast_results(input_file, query_lengths, reference_lengths, output_file, database_name):

    with open(input_file, 'r') as file:
        lines = file.readlines()

    # Filter out lines that contain '#'
    filtered_lines = [line for line in lines if '#' not in line]
    file_path = './Results/' + database_name + '/blast_processed.out'
    with open(file_path, 'w') as file:
        file.writelines(filtered_lines)

    # Read the BLAST result file
    file_path = './Results/' + database_name + '/blast_processed.out'
    df = pd.read_csv(file_path, sep='\t', header=None)

    df.columns = ['query_id', 'reference_id', 'identity', 'alignment_length', 'mismatches', 'gap opens', 'q start', 'q end', 's start', 's end', 'e value', 'bit score']
    
    # 比对长度计算（区间合并）
    result = df.groupby(['query_id', 'reference_id']).apply(
        lambda x: pd.Series({
            'intervals': [(row['s start'], row['s end']) for _, row in x.iterrows()]
        })
    ).reset_index()
    
    result['interval_length'] = result['intervals'].apply(lambda x: merge_intervals(x))
    result = result.drop(columns=['intervals'])

    # Take the first four columns 
    df= df.iloc[:,:4 ]
    df.columns = ['query_id', 'reference_id', 'similarity', 'alignment_length']
    
    # The sum of the weighted similarity level and comparison length is calculated, grouped by query_id and reference_id
    grouped = df.groupby(['query_id', 'reference_id']).apply(lambda x: pd.Series({
        'similarity': (x['similarity'] * x['alignment_length']).sum() / x['alignment_length'].sum()
    })).reset_index()
    
    # Read length information for query and reference
    query_lengths.columns = ['query_id', 'query_length']
    
    reference_lengths.columns = ['reference_id', 'reference_length']
    
    # Merge length information to the grouped data frame
    grouped = pd.merge(grouped, query_lengths, on='query_id', how='left')
    grouped = pd.merge(grouped, reference_lengths, on='reference_id', how='left')
    grouped = pd.merge(grouped, result, on=['query_id', 'reference_id'], how='left')

    grouped['query_length'] = grouped['query_length'].astype(float)
    grouped['reference_length'] = grouped['reference_length'].astype(float)
    
    grouped['alignment_length/query_length'] = grouped['interval_length'] / grouped['query_length']
    grouped['alignment_length/reference_length'] = grouped['interval_length'] / grouped['reference_length']

    # Save the result to a new file
    grouped.to_csv(output_file, sep='\t', header=False, index=False)

# Integrate tmvec data and blast data
def tmvec_corresponds_to_blast(tmvec_files, Blast_merge_files, query_lengths, reference_lengths, database_name):
    
    df1 = pd.read_csv(tmvec_files, sep='\t')

    df2 = pd.read_csv(Blast_merge_files, sep='\t', header=None)

    df2.columns = ['query_id', 'reference_id', 'similarity', 'query_length', 'reference_length', 'interval_length', 'alignment_length/query_length', 'alignment_length/reference_length']
    
    nan_df = pd.DataFrame({'similarity': 0, 'interval_length': 0, 'alignment_length/query_length' : 0, 'alignment_length/reference_length' : 0}, index=[1])

    df3 = pd.DataFrame()
    for i in range(len(df1)):
        query_df = df2[df2['query_id'].str.contains(df1['query_id'][i], case=True, na= False, regex=False)]
        query_database_df = query_df[query_df['reference_id'].str.contains(df1['database_id'][i], case=True, na= False, regex=False)]
        if len(query_database_df) > 0:
            df3 = pd.concat([df3, query_database_df[['similarity','interval_length', 'alignment_length/query_length', 'alignment_length/reference_length']][0:1]], ignore_index=True)
        else:
            df3 = pd.concat([df3, nan_df], ignore_index=True)
            
    tmvec_blast_df = pd.concat([df1, df3], axis = 1)


    query_lengths.columns = ['query_id', 'query_length']
    query_lengths = query_lengths.drop_duplicates('query_id') # 序列文件中有重复，id相同的序列
    tmvec_blast_df = pd.merge(tmvec_blast_df, query_lengths, on='query_id', how='left')
    # 得到reference长度文件
    reference_lengths.columns = ['database_id', 'reference_length']
    tmvec_blast_df = pd.merge(tmvec_blast_df, reference_lengths, on='database_id', how='left')
    output_file='./Results/' + database_name + '/tmvec_blast_df.tmvec'
    tmvec_blast_df.to_csv(output_file, sep='\t')

    # Optional
    # 筛选出符合要求的基因对 1为相似 0为不相似
    # 情况1：sequence(0) structure(1) length(1)
    # 情况2：sequence(0) structure(1) length(0)
    # 情况3：sequence(1) structure(1) length(0)
    # 情况4：sequence(1) structure(0) length(1)
    similarity_control = tmvec_blast_df['similarity'] < 50
    tm_score_control = tmvec_blast_df['tm-score'] > 0.5
    # query_seq_length_control = tmvec_blast_df['alignment_length/query_length'] > 0.5
    # seq_length_control = ((tmvec_blast_df['query_length'] / tmvec_blast_df['reference_length']) > 0.9) & ((tmvec_blast_df['query_length'] / tmvec_blast_df['reference_length']) < 1.1)

    # reference_seq_length_control = tmvec_blast_df['alignment_length/reference_length'] > 0.5
    filter_gene_df = tmvec_blast_df[similarity_control & tm_score_control]



    # Screen key gene pairs for nitrogen fixation (From: Celebrating 20 Years of Genetic Discoveries in Legume Nodulation and Symbiotic Nitrogen Fixation)
    key_node_genes_sym = pd.DataFrame({
        'gene_symbol' : ['MtIPD3L', 'MtDMI1(LjPOLLUX)', 'LjCASTOR', 'LjSYMRK', 'MtCCaMK(LOF)', 'MtHMGR1', 'LjLNP', 'LjNFR1_MtLYK3', 'LjNFR5', 'MtNSP1', 'MtNSP2', 'LjnsRING', 'LjNUP85', 'PvRabA2']
        })
    filtered_key_nodule_df = filter_gene_df[filter_gene_df['query_id'].apply(lambda x: any(gene in str(x) for gene in key_node_genes_sym['gene_symbol']))]

    filtered_key_nodule_df['query_gene_name'] = filtered_key_nodule_df['query_id'].str.split('__').str[0]
    output_file='./Results/' + database_name + '/filtered_key_nodule_df.tmvec'
    filtered_key_nodule_df.to_csv(output_file, sep='\t')

    return filtered_key_nodule_df