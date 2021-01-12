rmdir -rf ./obj
mkdir ./obj

g++ -c -O3 others.cpp -o ./obj/others.o
g++ -c -O3 mwis.cpp -o ./obj/mwis.o
g++ mwis.o others.o -o mwis