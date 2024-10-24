# Benchmarking MSA

## Aligners:
1. FAMSA (bash and python)
2. Kalign3
3. MAFFT-PartTree (Commented out; memory consuming)
4. ClustalO (Commented out; time consuming)

## Dataset:
Download extHomFam-v2 from zenodo [link](https://zenodo.org/records/6524237) to project directory.

```bash
wget https://zenodo.org/records/6524237/files/extHomFam-v2.zip
unzip -qq extHomFam-v2.zip
rm extHomFam-v2.zip
```

# Example Usage
Using synthetic data from extHomFam-v2 medium:
```bash
THREADS=128
python3 ./benchmark-MSA/benchmark.py \
        --threads ${THREADS}
```

Using all data from extHomFam-v2 (small, medium, large, xlarge):
```bash
THREADS=128
python3 ./benchmark-MSA/benchmark.py \
        --threads ${THREADS} \
        --whole-extHomFam-v2
```
