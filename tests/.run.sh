touch tests//.timeout
CMD="valgrind --leak-check=full /home/alvaro/UGRMPBase/NetBeansProjects/MiKmer_4/dist/Debug/GNU-Linux/mikmer_4  ../Genomes/human1_k1.prf ../Genomes/human1.prf ../Genomes/human2.prf 1> tests//.out11 2>&1"
eval $CMD
rm tests//.timeout
