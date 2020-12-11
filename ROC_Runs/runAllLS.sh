
cp ./updateROC.R adult/ 
cp ./updateROC.R dauer/
cp ./updateROC.R L1/
cp ./updateROC.R L2/
cp ./updateROC.R L3/
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
cd  "/data1/compbio/zlu21/AcrossTissue/RunResults/L1"
./run.sh &
cd  "/data1/compbio/zlu21/AcrossTissue/RunResults/L2"
./run.sh &
cd  "/data1/compbio/zlu21/AcrossTissue/RunResults/L3"
./run.sh &
cd "/data1/compbio/zlu21/AcrossTissue/RunResults/emb"
./run.sh &



