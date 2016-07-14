CC='g++'
LIB='-lgsl -lgslcblas'
OPT='-O2'
TARGET='main'
SRC='main.cpp'

$CC $OPT -o $TARGET $SRC $LIB
