
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



cd "/data1/compbio/zlu21/AcrossTissue/RunResults/adult"
./run.sh &
cd "/data1/compbio/zlu21/AcrossTissue/RunResults/dauer"
./run.sh &
cd  "/data1/compbio/zlu21/AcrossTissue/RunResults/L3L4"
./run.sh &
cd "/data1/compbio/zlu21/AcrossTissue/RunResults/larvae"
./run.sh &
cd "/data1/compbio/zlu21/AcrossTissue/RunResults/emb"
./run.sh &



