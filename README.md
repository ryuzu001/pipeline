Pipeline project for UCR's CS130 course. 

after compilation, run minigl or minigl-nogl with the following flags:

Usage: ./minigl-nogl -i <input-file> [ -s <solution-file> ] [ -o <stats-file> ] [ -c ] [ -p ] [ -g ]
    <input-file>      File with commands to run
    <solution-file>   File with solution to compare with
    <stats-file>      Dump statistics to this file rather than stdout
    -p                Dump results to png files
    -g                Use OpenGL
    -c                Compare only; if using OpenGL, don't leave window open
    -n                Run solution multiple times to get more reliable timing results

For example, to run an individual test without openGL, use 

    ./minigl-nogl -i tests/xx.txt -s tests/xx.png -p

where xx is the number of the specific test to be used. 

    ./minigl-nogl -i tests/00.txt -s tests/00.png -p

You can compile with scons if installed or make otherwise.

To run the grading script, use 

    ./grading-script.sh ./tests