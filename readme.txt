CSCI 716: Assignment 3 (Trapezoidal Map) -- Yuying Mao & Owen Shriver

trapezoidalmap.py is the python program implementing the Random Incremental Algorithm for solving trapezoidal mapping problem. It reads the bounding box and line segments from an inputfile, performs the algorithm, and outputs the corresponding matrix along with a plot of the map.

usage:
    python trapezoidalmap.py inputfile

The program creates two files: matrix.txt, which contains the adjacency matrix, and plot.png, which contains an image of the trapezoidal map.

Then, the user may input a point in the form "(x, y)" (with the parentheses and comma), and the program will print to standard output the path to that point in the map.

To make the results visualized, a library call 'matplotlib' is imported in this program
For more information: https://matplotlib.org/index.html

to install matplotlib:
    python -m pip install -U matplotlib