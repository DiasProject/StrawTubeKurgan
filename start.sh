rm -rf CMake*
cmake ../B1
make -j4
./exampleB1 begin_mu.g4 > test2.txt
# delete between lines $a and $b inclusive
more test2.txt | awk '0 <= NR && NR <= 734 {next} {print}' > test.txt
rm test2.txt
#more test1.txt | awk 'FNR==0{print $0; next} {print FNR,$0}' OFS=" " > test.txt
#rm test1.txt
