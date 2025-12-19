g++ -std=c++11 `root-config --cflags` -Ilib -Llib -c sum.cc
g++ -Ilib -Llib -std=legacy -o sum.exe sum.o `root-config --libs`

g++ -std=c++11 `root-config --cflags` -Ilib -Llib -c sum_mc.cc
g++ -Ilib -Llib -std=legacy -o sum_mc.exe sum_mc.o `root-config --libs`
