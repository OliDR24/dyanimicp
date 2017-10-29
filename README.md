# Dynamic-Programming

You are required to do your assignment of dynamic programming within this repository and your last commit before the deadline will be counted for grading, though you could keep pushing your code after the deadline. 

## Requirements of the repository

`exchange_matrices` is a folder containing all the exchange matrices used in the script `align.py`. **Please don't modify or relocate it**.

`align.py` is the skeleton script where you need to add 'flesh' to make it work. **You are not allowed to rename or relocate the script.** It is possible to have multiple scripts, but `align.py` has to be the main script.

`student.id` is a text file to store your personal information to help us identify whose repository belongs to whom. In this file, you have found the student information of Cico Zhang, you are required to change it to yours. **The first line is your full name; second line, student ID; third line, VUnet ID.**

## Requirements of the script `align.py`

### Run this script

Please guarantee that the command `python3 align.py -f [input file] -e [exchange matrices] -l (or -g or -s) -p [penalty] -o [output] -m [score matrix]` always work no matter how you customise the script `align.py`. (This repository is the home directory when running the script). It also means that your code has to be at least compatible with Python 3.

| parameter | Explanation|
|:---------:|:----------:|
|`-f`| designates input file |
|`-o`| designates the path and name of the file to store the alignments|
|`-m`| designates the path and name of the file to store the score matrix|
|`-l` or `-g` or `-s`| designates the algorithm: local, global or semi-global|
|`-e`| designates the exchange matrix: identity, pam250, blosum62|
|`-p`| designates the penalty, 2 is the default value|
|`-v`| prints the output to the screen; this will make your testing easy|

For example,
> python3 align.py -f test.fasta -e pam250 -l -p 3 -o test.align -m test.scorematrix

## The output format

There are two major outputs for each execution of the script: the alignment and the score matrix. 

The alignment(s) and score matrix are stored in a list of lists in Python, in which space, letters, `|` and numbers are the possible elements. 

In the skeleton script, we have provided the function `print_matrix_to_file` that writes the alignment and score matrix into files (for grading), we have also provided the function `print_matrix_on_screen` that prints the alignment and score matrix on the screen to facilitate your testing.
