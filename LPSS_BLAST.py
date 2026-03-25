from utils import *
import pandas as pd
import os
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--query_gene", type=str, required=True)
    parser.add_argument("--database_gene", type=str, required=True)
    args = parser.parse_args()
    print(args)
database_name = os.path.splitext(os.path.basename(args.database_gene))[0]

# Blast
print('Start the Blast process.')
makeblastdb_cmd = 'makeblastdb -in ' + args.database_gene + ' -dbtype prot -out ./DB/' + database_name + '/blast_database'
if os.path.exists('./Results/' + database_name):
    print('The results folder already exists.')
else:
    os.makedirs('./Results/' + database_name)

blastp_cmd = 'blastp -db ./DB/' + database_name + '/blast_database -query ' + args.query_gene + ' -out ./Results/' + database_name + '/blast.out -outfmt 7 -evalue 1e-5'
if os.path.exists('./Results/' + database_name + '/blast.out'):
    print('The blast result already exists.')
elif os.path.exists('./DB/' + database_name + '/blast_database.pdb'):
    print('The blast database already exists.')
    os.system(blastp_cmd)
else:
    os.system(makeblastdb_cmd)
    os.system(blastp_cmd)


# query length
query_lengths = get_length_file(args.query_gene)

# reference length
reference_lengths = get_length_file(args.database_gene)

input_file = './Results/' + database_name + '/blast.out'
output_file = './Results/' + database_name + '/blast.merge'
process_blast_results(input_file, query_lengths, reference_lengths, output_file, database_name)

# Tm-vec
print('Start the Tm-vec process.')
maketmvecdb_cmd = 'tmvec-build-database --input-fasta ' + args.database_gene + ' --tm-vec-model ./model/tm_vec_swiss_model.ckpt --tm-vec-config-path ./model/tm_vec_swiss_model_params.json --device \'gpu\' --output ./DB/' + database_name + '/tmvec_database'
tmvec_search_cmd = 'tmvec-search --query ' + args.query_gene + ' --tm-vec-model ./model/tm_vec_swiss_model.ckpt --tm-vec-config ./model/tm_vec_swiss_model_params.json --database ./DB/' + database_name + '/tmvec_database/db.npy --metadata ./DB/' + database_name + '/tmvec_database/meta.npy --database-fasta ' + args.query_gene + ' --device \'gpu\' --output-format tabular --output ./Results/' + database_name + '/tabular.txt --output-embeddings ./Results/' + database_name + '/test.npy --k-nearest-neighbors 50'
if os.path.exists('./Results/' + database_name + '/tabular.txt'):
    print('The tmvec result already exists.')
elif os.path.exists('./DB/' + database_name + '/tmvec_database'):
    print('The tmvec database already exists.')
    os.system(tmvec_search_cmd)
else:
    os.system(maketmvecdb_cmd)
    os.system(tmvec_search_cmd)

# Corresponding and filtering
print('Start the corresponding process.')
tmvec_corresponds_to_blast_filter_df = tmvec_corresponds_to_blast('./Results/' + database_name + '/tabular.txt', './Results/' + database_name + '/blast.merge', query_lengths, reference_lengths, database_name)
print(tmvec_corresponds_to_blast_filter_df)
print('Finished.')