python gprof2dot.py -e 0.0 profile.out > prog.dot
dot -Tpng -o profile_euler.png prog.dot
display profile_euler.png
