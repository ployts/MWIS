mkdir ./obj

g++ -c -O3 ./src/others.cpp -o ./obj/others.o
g++ -c -O3 ./src/mwis.cpp -o ./obj/mwis.o
g++ -c -O3 ./src/main.cpp -o ./obj/main.o
g++ ./obj/*.o -o mwis

rm -rf ./obj    