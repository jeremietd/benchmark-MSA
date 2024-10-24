import argparse
from aligners import famsa, famsa_medoid, clustalo, mafft_parttree, kalign3, famsa_python, famsa_medoid_python
from util import prepare_extHomFam_v2, parse_fasta, create_synthetic_dataset, dict_to_dataframe
import os
import pandas as pd
import tqdm
import tempfile as tmp

# Define arguments
parser = argparse.ArgumentParser()
parser.add_argument('--threads', type=int, help='Number of threads', required=True)
parser.add_argument('--whole-extHomFam-v2', action='store_true', help='Use the whole extHomFam-v2 dataset')
parser.add_argument('--no-python', action='store_true', help='Do not use the Python implementation of FAMSA and FAMSA-Medoid')
args = parser.parse_args()

result_dict = {}
args.threads = str(args.threads)

# Prepare extHomFam-v2 dataset
extHomFam_v2 = prepare_extHomFam_v2(all=args.whole_extHomFam_v2)

# Take the xlarge dataset and randomly take without replacing 100,000 (small), 250,000 (medium), 500,000 (large) sequences
mini_extHomFam_v2 = create_synthetic_dataset(extHomFam_v2)

dataset_for_use = extHomFam_v2 if args.whole_extHomFam_v2 else mini_extHomFam_v2

#Prepare save path
print(f"Skipping Python implementation of FAMSA and FAMSA-Medoid: {args.no_python}") if args.no_python else print("Using Python implementation of FAMSA and FAMSA-Medoid")
with tmp.TemporaryDirectory() as tmpdirname:
        print(f"Created temporary directory: {tmpdirname}")
        save_path = tmpdirname
        # Benchmarking
        pbar = tqdm.tqdm(total=len(dataset_for_use.keys()))
        for sizes in list(dataset_for_use.keys()):
                file_name = dataset_for_use[sizes]
                clean_file_name = file_name.split("/")[-1].split(".")[0]
                print(f"File name: {file_name}")

                print(f"=================ALIGNING {sizes.upper()}=================")
                famsa(input_file=f"{file_name}", 
                        output_file=f"{save_path}/{clean_file_name}-FAMSA.fasta", 
                        threads=args.threads,
                        dataset_size=sizes,
                        result_dict=result_dict)
                print(f"Complete FAMSA")
                
                famsa_medoid(input_file=f"{file_name}",
                        output_file=f"{save_path}/{clean_file_name}-FAMSA-medoid.fasta",
                        threads=args.threads,
                        dataset_size=sizes,
                        result_dict=result_dict)
                print(f"Complete FAMSA-Medoid")

                if not args.no_python:
                        famsa_python(input_file=f"{file_name}",
                                output_file=f"{save_path}/{clean_file_name}-FAMSA-python.fasta",
                                threads=1,
                                dataset_size=sizes,
                                result_dict=result_dict)
                        print(f"Complete FAMSA-Python")
                
                        famsa_medoid_python(input_file=f"{file_name}",
                                output_file=f"{save_path}/{clean_file_name}-FAMSA-medoid-python.fasta",
                                threads=1,
                                dataset_size=sizes,
                                result_dict=result_dict)
                        print(f"Complete FAMSA-Medoid-Python")
                
                # #CLUSTALO # takes long time
                # clustalo(input_file=f"{file_name}",
                #         output_file=f"{save_path}/{clean_file_name}-CLUSTALO.fasta",
                #         threads=args.threads,
                #         dataset_size=sizes,
                #         result_dict=result_dict)
                # print(f"Complete CLUSTALO")
                
                # # MAFFT-PartTree # uses too much memory
                # mafft_parttree(input_file=f"{file_name}",
                #         output_file=f"{save_path}/{clean_file_name}-MAFFT.fasta",
                #         threads=args.threads,
                #         dataset_size=sizes,
                #         result_dict=result_dict)
                # print(f"Complete MAFFT-PartTree")
                
                kalign3(input_file=f"{file_name}",
                        output_file=f"{save_path}/{clean_file_name}-KALIGN3.fasta",
                        threads=args.threads,
                        dataset_size=sizes,
                        result_dict=result_dict)
                print(f"Complete KALIGN3")
                
                pbar.update(1)
                print(f"Aligned {sizes.upper()}")
        pbar.close()

df = dict_to_dataframe(result_dict)
folder_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "MSAresults")
os.makedirs(folder_path, exist_ok=True)

if args.whole_extHomFam_v2: df.to_csv(f"{folder_path}/MSA-results_{args.threads}_whole.csv")
else: df.to_csv(f"{folder_path}/MSA-results_{args.threads}_synthetic.csv")
