
cp ./updateROC.R adult/ 
cp ./updateROC.R dauer/
cp ./updateROC.R L3L4/
cp ./updateROC.R larvae/
cp ./updateROC.R emb/

cp ./run.sh adult/
cp ./run.sh dauer/
cp ./run.sh L3L4/
cp ./run.sh larvae/
cp ./run.sh emb/



adult/run.sh &
dauer/run.sh &
L3L4/run.sh &
larvae/run.sh &
emb/run.sh &



