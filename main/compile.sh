rm *.o event_eval.exe 
g++ -std=c++11 `root-config --cflags` -Iinclude -I/eic/u/cvhulse/epic/lhapdf/include/ -Ilib -Llib -c event_eval.cc
gfortran -std=legacy -I/usr/include -I/eic/u/cvhulse/epic/lhapdf/include/ -lgfortran -L/usr/lib/ -L/eic/u/cvhulse/epic/lhapdf/lib/ -lHAPDF -c DSSV.f vegas.f fDSS.f fuu.f deltaq.f
g++ -I./include -I/eic/u/cvhulse/epic/lhapdf/include -I/eic/u/cvhulse/epic/lhapdf/share -Ilib -Llib -std=legacy -o event_eval.exe event_eval.o DSSV.o vegas.o fDSS.o fuu.o deltaq.o `root-config --libs` -lgfortran -L/direct/eic+u/cvhulse/ecce/lhapdf/lib/ -lLHAPDF
