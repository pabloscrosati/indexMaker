# Index Maker for Finding and Generating Protein-Solute Contacts in GROMACS MD Simulations
**Written by: Pablo M. Scrosati**
**Version 0.9.1 beta**

"beta" designation will be kept until the program has been externally validated.

## Usage
`main.py [-h] -f  -d  -o  [--csv] [--st] [--coarse-cutoff] [--short-cutoff]`

```
  -h , --help            show this help message and exit
  -f , --file-pattern   filename pattern for input files
  -d , --directory      directory containing files to import
  -o , --output         output file
  --csv                 output CSV file alongside tab-delimited file
  --st                  run program in single-threaded mode
  --coarse-cutoff       define coarse cutoff distance
  --short-cutoff        define short cutoff distance
```

## Folder Structure for Usage
```
parent_directory
├─── main.py
├────── <directory>
│       ├─── <file_pattern>_1.gro
│       ├─── <file_pattern>_2.gro
│       ├─── <file_pattern>_3.gro
│       ├─── ···
│       └─── <file_pattern>_n.gro
└────── output (will be created if not present)
        ├─── <output>.txt
        ├─── <output>.ndx
        └─── [<output>.csv]
```
*Note: Individual file reports are provided within `<directory>`.*


<!-- │ ├ ─ └ -->