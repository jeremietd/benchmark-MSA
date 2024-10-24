import glob
import os
import random
import pandas as pd

folder_path = os.path.dirname(os.path.realpath(__file__))

def concat_fasta(input_file, output_file):
  with open(output_file, 'w') as outfile:
    for file in glob.glob(f"{input_file}/*"):
      with open(file) as infile:
        outfile.write(infile.read())

def prepare_extHomFam_v2(all=False):
    check_list = ["extHomFam-v2-small.fasta", "extHomFam-v2-medium.fasta", "extHomFam-v2-large.fasta", "extHomFam-v2-xlarge.fasta"]
    # https://zenodo.org/records/6524237
    if all and all([os.path.exists(f"{folder_path}/{check_list[i]}") for i in range(4)]):
        return {
            "small": f"{folder_path}/extHomFam-v2-small.fasta",
            "medium": f"{folder_path}/extHomFam-v2-medium.fasta",
            "large": f"{folder_path}/extHomFam-v2-large.fasta",
            "xlarge": f"{folder_path}/extHomFam-v2-xlarge.fasta"
            }
    elif os.path.exists(f"{folder_path}/{check_list[1]}"):
        return {
            "medium": f"{folder_path}/extHomFam-v2-medium.fasta"
            }

    else:
        if all:
            concat_fasta(f"{folder_path}/extHomFam-v2/small", f"{folder_path}/extHomFam-v2-small.fasta")
            concat_fasta(f"{folder_path}/extHomFam-v2/medium", f"{folder_path}/extHomFam-v2-medium.fasta")
            concat_fasta(f"{folder_path}/extHomFam-v2/large", f"{folder_path}/extHomFam-v2-large.fasta")
            concat_fasta(f"{folder_path}/extHomFam-v2/xlarge", f"{folder_path}/extHomFam-v2-xlarge.fasta")

            return {
            "small": f"{folder_path}/extHomFam-v2-small.fasta",
            "medium": f"{folder_path}/extHomFam-v2-medium.fasta",
            "large": f"{folder_path}/extHomFam-v2-large.fasta",
            "xlarge": f"{folder_path}/extHomFam-v2-xlarge.fasta"
            }

        else:
            concat_fasta(f"{folder_path}/extHomFam-v2/medium", f"{folder_path}/extHomFam-v2-medium.fasta")

            return {"medium": f"{folder_path}/extHomFam-v2-medium.fasta"}
        

def create_synthetic_dataset(extHomFam_v2):
  output = {
    "xsmall": f"{folder_path}/mini-extHomFam-v2-xsmall.fasta",
    "small": f"{folder_path}/mini-extHomFam-v2-small.fasta", 
    "medium": f"{folder_path}/mini-extHomFam-v2-medium.fasta", 
    "large": f"{folder_path}/mini-extHomFam-v2-large.fasta"}

  # check if all output files exist
  if all([os.path.exists(output[key]) for key in output.keys()]):
    return output
  else:
    sequences = []
    for name, seq in zip(*parse_fasta(extHomFam_v2["medium"], return_names=True, clean=None, full_name=False)):
      sequences.append((name,seq))
    random.shuffle(sequences)

    base_size = 50000
    mini_extHomFam_v2_xsmall = sequences[:base_size] # 50,000
    mini_extHomFam_v2_small = sequences[base_size:base_size+100000] # 100,000
    mini_extHomFam_v2_medium = sequences[base_size+100000:base_size+350000] # 250,000
    mini_extHomFam_v2_large = sequences[base_size+350000:base_size+850000] # 500,000

    with open(f"{folder_path}/mini-extHomFam-v2-xsmall.fasta", "w") as xsmall, open(f"{folder_path}/mini-extHomFam-v2-small.fasta", "w") as small, open(f"{folder_path}/mini-extHomFam-v2-medium.fasta", "w") as medium, open(f"{folder_path}/mini-extHomFam-v2-large.fasta", "w") as large:
            print(f"Creating XSMALL dataset with {len(mini_extHomFam_v2_xsmall)} sequences")
            for name, seq in mini_extHomFam_v2_xsmall:
                    print(f">{name}\n{seq}", file=xsmall)
            
            print(f"Creating SMALL dataset with {len(mini_extHomFam_v2_small)} sequences")
            for name, seq in mini_extHomFam_v2_small:
                    print(f">{name}\n{seq}", file=small)

            print(f"Creating MEDIUM dataset with {len(mini_extHomFam_v2_medium)} sequences")
            for name, seq in mini_extHomFam_v2_medium:
                    print(f">{name}\n{seq}", file=medium)
            
            print(f"Creating LARGE dataset with {len(mini_extHomFam_v2_large)} sequences")
            for name, seq in mini_extHomFam_v2_large:
                    print(f">{name}\n{seq}", file=large)
    
    return output
                
def save_results(result_dict, aligner, current, peak, diff, time, dataset_size, threads):
    if aligner not in result_dict:
        result_dict[aligner] = {}

    result_dict[aligner][dataset_size] = {
        "current": current,
        "peak": peak,
        "usage": diff,
        "time": time,
        "threads": threads
    }

    return result_dict


def dict_to_dataframe(result_dict):
    data = []

    for aligner, dataset_sizes in result_dict.items():
        for dataset_size, values in dataset_sizes.items():
            row = [aligner, dataset_size] + list(values.values())
            data.append(row)
    
    columns = ['Aligner', 'Dataset Size', 'Current', 'Peak', 'Usage', 'Time', 'Threads']
    
    df = pd.DataFrame(data, columns=columns)
    
    return df
    
def _open_if_is_name(filename_or_handle, mode="r"):
    """
        if a file handle is passed, return the file handle
        if a Path object or path string is passed, open and return a file handle to the file.

        returns:
            file_handle, input_type ("name" | "handle")
    """
    out = filename_or_handle
    input_type = "handle"
    try:
        out = open(filename_or_handle, mode)
        input_type = "name"
    except TypeError:
        pass
    except Exception as e:
        raise(e)

    return (out, input_type)

def parse_fasta(filename, return_names=False, clean=None, full_name=False): 
    """
        adapted from: https://bitbucket.org/seanrjohnson/srj_chembiolib/src/master/parsers.py
        

        input:
            filename: the name of a fasta file or a filehandle to a fasta file.
            return_names: if True then return two lists: (names, sequences), otherwise just return list of sequences
            clean: {None, 'upper', 'delete', 'unalign'}
                    if 'delete' then delete all lowercase "." and "*" characters. This is usually if the input is an a2m file and you don't want to preserve the original length.
                    if 'upper' then delete "*" characters, convert lowercase to upper case, and "." to "-"
                    if 'unalign' then convert to upper, delete ".", "*", "-"
            full_name: if True, then returns the entire name. By default only the part before the first whitespace is returned.

        output: sequences or (names, sequences)
    """
    
    prev_len = 0
    prev_name = None
    prev_seq = ""
    out_seqs = list()
    out_names = list()
    (input_handle, input_type) = _open_if_is_name(filename)

    for line in input_handle:
        line = line.strip()
        if len(line) == 0:
            continue
        if line[0] == ">":
            if full_name:
                name = line[1:]
            else:
                line = line.replace('|', '/')
                parts = line.split('/', 1)
                name = parts[0][1:]
            out_names.append(name)
            if (prev_name is not None):
                out_seqs.append(prev_seq)
            prev_len = 0
            prev_name = name
            prev_seq = ""
        else:
            prev_len += len(line)
            prev_seq += line
    if (prev_name != None):
        out_seqs.append(prev_seq)

    if input_type == "name":
        input_handle.close()
    
    if clean == 'delete':
        # uses code from: https://github.com/facebookresearch/esm/blob/master/examples/contact_prediction.ipynb
        deletekeys = dict.fromkeys(string.ascii_lowercase)
        deletekeys["."] = None
        deletekeys["*"] = None
        translation = str.maketrans(deletekeys)
        remove_insertions = lambda x: x.translate(translation)

        for i in range(len(out_seqs)):
            out_seqs[i] = remove_insertions(out_seqs[i])
    
    elif clean == 'upper':
        deletekeys = {'*': None, ".": "-"}
        translation = str.maketrans(deletekeys)
        remove_insertions = lambda x: x.translate(translation)

        for i in range(len(out_seqs)):
            out_seqs[i] = remove_insertions(out_seqs[i].upper())
    elif clean == 'unalign':
        anyX = str(random.choice(ESM_ALLOWED_AMINO_ACIDS))
        anyB = str(random.choice("ND"))
        anyZ = str(random.choice("EQ"))
        deletekeys = {'*': None, ".": None, "-": None, "X": anyX, "B": anyB, "Z": anyZ}
        
        translation = str.maketrans(deletekeys)
        remove_insertions = lambda x: x.translate(translation)
        
        for i in range(len(out_seqs)):
            out_seqs[i] = remove_insertions(out_seqs[i].upper())
    elif clean is not None:
        raise ValueError(f"unrecognized input for clean parameter: {clean}")

    if return_names:
        return out_names, out_seqs
    else:
        return out_seqs