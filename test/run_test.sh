make # ensure heptagon build
mkdir -p test/out_trees/

simname=sim_close$1

printf "======================\n\n\e[0;32mRunning HEPTAGON\e[0m\n"
time ./heptagon --in test/$simname.dist.phylip --skeleton-size 30000 > test/out_trees/heptagon.txt
python3 MAPLE.py --inputTree test/$simname.treefile --inputRFtrees test/out_trees/heptagon.txt --output test/out_trees/ > /dev/null
echo $'\nNormalized Robinson Foulds Distance:'
cat test/out_trees/_RFdistances.txt | grep -oE '\s0\.[0-9]+' | head -n1
echo $'\nNormalized Quartet Distance:'
./quartet_dist -v test/$simname.treefile test/out_trees/heptagon.txt | grep -oE '\s0\.[0-9]+' | head -n1
rm test/out_trees/heptagon.txt
rm test/out_trees/_RFdistances.txt

if [ -f quicktree ] ; then
printf "======================\n\n\e[0;32mRunning HEPTAGON\e[0m\n"
time ./quicktree -in m test/$simname.dist.phylip | tr -d '\n' > test/out_trees/quicktree.txt
python3 MAPLE.py --inputTree test/$simname.treefile --inputRFtrees test/out_trees/quicktree.txt --output test/out_trees/ > /dev/null
echo $'\nNormalized Robinson Foulds Distance:'
cat test/out_trees/_RFdistances.txt | grep -oE '\s0\.[0-9]+' | head -n1
echo $'\nNormalized Quartet Distance:'
./quartet_dist -v test/$simname.treefile test/out_trees/quicktree.txt | grep -oE '\s0\.[0-9]+' | head -n1
rm test/out_trees/quicktree.txt
rm test/out_trees/_RFdistances.txt
fi;

if [ -f dipper ] ; then
printf "======================\n\n\e[0;32mRunning HEPTAGON\e[0m\n"
time ./dipper -i d -o t -I test/$simname.dist.phylip -O test/out_trees/dipper.txt
python3 MAPLE.py --inputTree test/$simname.treefile --inputRFtrees test/out_trees/dipper.txt --output test/out_trees/ > /dev/null
echo $'\nNormalized Robinson Foulds Distance:'
cat test/out_trees/_RFdistances.txt | grep -oE '\s0\.[0-9]+' | head -n1
echo $'\nNormalized Quartet Distance:'
./quartet_dist -v test/$simname.treefile test/out_trees/dipper.txt | grep -oE '\s0\.[0-9]+' | head -n1
rm test/out_trees/dipper.txt
rm test/out_trees/_RFdistances.txt   
fi;
