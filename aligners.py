import subprocess
import time
import tracemalloc
from util import save_results
from pyfamsa import Aligner, Sequence
from util import parse_fasta

def famsa_python(input_file, output_file, threads, dataset_size, result_dict):

    sequences = []
    for name, seq in zip(*parse_fasta(input_file, return_names=True, clean=None, full_name=False)):
        sequences.append(Sequence(name.encode(), seq.encode()))
    
    # start computing time and memory
    tracemalloc.start()
    start = time.time()

    aligner = Aligner(guide_tree="sl")
    msa = aligner.align(sequences)

    with open(output_file, "w") as f:
        for seq in msa:
            print(f">{seq.id.decode()}\n{seq.sequence.decode()}", file=f)

    end = time.time()
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    # Metrics
    time_taken = (end - start) / 60
    current_mb = current / 10**3
    peak_mb = peak / 10**3
    usage_mb = peak_mb - current_mb
    # print(f"Current memory usage is {current / 10**3}KB; Peak was {peak / 10**3}KB; Diff = {(peak - current) / 10**3}KB")
    # print(f"Total time taken: {total_time}")

    save_results(result_dict, "famsa-python", current_mb, peak_mb, usage_mb, time_taken, dataset_size, threads)

def famsa_medoid_python(input_file, output_file, threads, dataset_size, result_dict):

    sequences = []
    for name, seq in zip(*parse_fasta(input_file, return_names=True, clean=None, full_name=False)):
        sequences.append(Sequence(name.encode(), seq.encode()))
    
    # start computing time and memory
    tracemalloc.start()
    start = time.time()

    aligner = Aligner(guide_tree="sl", tree_heuristic="medoid")
    msa = aligner.align(sequences)

    with open(output_file, "w") as f:
        for seq in msa:
            print(f">{seq.id.decode()}\n{seq.sequence.decode()}", file=f)

    end = time.time()
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    # Metrics
    time_taken = (end - start) / 60
    current_mb = current / 10**3
    peak_mb = peak / 10**3
    usage_mb = peak_mb - current_mb
    # print(f"Current memory usage is {current / 10**3}KB; Peak was {peak / 10**3}KB; Diff = {(peak - current) / 10**3}KB")
    # print(f"Total time taken: {total_time}")

    save_results(result_dict, "famsa-medoid-python", current_mb, peak_mb, usage_mb, time_taken, dataset_size, threads)


def famsa(input_file, output_file, threads, dataset_size, result_dict):
    tracemalloc.start()
    start = time.time()

    try:
        proc = subprocess.run(["famsa", 
                                "-gz",
                                "-t", threads, 
                                input_file,
                                output_file],
                                check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                                )
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode('utf-8'))
        print(e.stdout.decode('utf-8'))
        raise e

    end = time.time()
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    # Metrics
    time_taken = (end - start) / 60
    current_mb = current / 10**3
    peak_mb = peak / 10**3
    usage_mb = peak_mb - current_mb
    # print(f"Current memory usage is {current / 10**3}KB; Peak was {peak / 10**3}KB; Diff = {(peak - current) / 10**3}KB")
    # print(f"Total time taken: {total_time}")

    save_results(result_dict, "famsa", current_mb, peak_mb, usage_mb, time_taken, dataset_size, threads)

def famsa_medoid(input_file, output_file, threads, dataset_size, result_dict):
    tracemalloc.start()
    start = time.time()

    try:
        proc = subprocess.run(["famsa", 
                                "-medoidtree",
                                "-gz",
                                "-t", threads, 
                                input_file,
                                output_file],
                                check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                                )
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode('utf-8'))
        print(e.stdout.decode('utf-8'))
        raise e

    end = time.time()
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    # Metrics
    time_taken = (end - start) / 60
    current_mb = current / 10**3
    peak_mb = peak / 10**3
    usage_mb = peak_mb - current_mb
    # print(f"Current memory usage is {current / 10**3}KB; Peak was {peak / 10**3}KB; Diff = {(peak - current) / 10**3}KB")
    # print(f"Total time taken: {total_time}")

    save_results(result_dict, "famsa-medoid", current_mb, peak_mb, usage_mb, time_taken, dataset_size, threads)

def clustalo(input_file, output_file, threads, dataset_size, result_dict):
    tracemalloc.start()
    start = time.time()

    try:
        proc = subprocess.run(["clustalo", 
                                "--threads", threads, 
                                "-i", input_file,
                                "-o", output_file,
                                "--force"],
                                check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                                )
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode('utf-8'))
        print(e.stdout.decode('utf-8'))
        raise e

    end = time.time()
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    # Metrics
    time_taken = (end - start) / 60
    current_mb = current / 10**3
    peak_mb = peak / 10**3
    usage_mb = peak_mb - current_mb
    # print(f"Current memory usage is {current / 10**3}KB; Peak was {peak / 10**3}KB; Diff = {(peak - current) / 10**3}KB")
    # print(f"Total time taken: {total_time}")

    save_results(result_dict, "clustalo", current_mb, peak_mb, usage_mb, time_taken, dataset_size, threads)

def mafft_parttree(input_file, output_file, threads, dataset_size, result_dict):
    tracemalloc.start()
    start = time.time()
    try:
        proc = subprocess.run(f"mafft --anysymbol --quiet --parttree --thread {threads} {input_file} > {output_file}",
                                shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode('utf-8'))
        print(e.stdout.decode('utf-8'))
        raise e

    end = time.time()
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    # Metrics
    time_taken = (end - start) / 60
    current_mb = current / 10**3
    peak_mb = peak / 10**3
    usage_mb = peak_mb - current_mb
    # print(f"Current memory usage is {current / 10**3}KB; Peak was {peak / 10**3}KB; Diff = {(peak - current) / 10**3}KB")
    # print(f"Total time taken: {total_time}")

    save_results(result_dict, "mafft-parttree", current_mb, peak_mb, usage_mb, time_taken, dataset_size, threads)

def kalign3(input_file, output_file, threads, dataset_size, result_dict):
    tracemalloc.start()
    start = time.time()

    try:
        proc = subprocess.run(["kalign", 
                                "--nthreads", threads, 
                                "-i", input_file,
                                "-o", output_file],
                                check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                                )
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode('utf-8'))
        print(e.stdout.decode('utf-8'))
        raise e

    end = time.time()
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    # Metrics
    time_taken = (end - start) / 60
    current_mb = current / 10**3
    peak_mb = peak / 10**3
    usage_mb = peak_mb - current_mb
    # print(f"Current memory usage is {current / 10**3}KB; Peak was {peak / 10**3}KB; Diff = {(peak - current) / 10**3}KB")
    # print(f"Total time taken: {total_time}")

    save_results(result_dict, "kalign3", current_mb, peak_mb, usage_mb, time_taken, dataset_size, threads)