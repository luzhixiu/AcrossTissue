
cp ./updateROC.R adult/ 
cp ./updateROC.R dauer/
cp ./updateROC.R L1/
cp ./updateROC.R L2L3/
cp ./updateROC.R emb/

cp ./run.sh adult/
cp ./run.sh dauer/
cp ./run.sh L1/
cp ./run.sh L2L3/
cp ./run.sh emb/



cd "/data1/compbio/zlu21/AcrossTissue/ROC_Runs/Run_6LS/adult"
./run.sh &
cd "/data1/compbio/zlu21/AcrossTissue/ROC_Runs/Run_6LS/dauer"
./run.sh &
cd  "/data1/compbio/zlu21/AcrossTissue/ROC_Runs/Run_6LS/L1"
./run.sh &
cd  "/data1/compbio/zlu21/AcrossTissue/ROC_Runs/Run_6LS/L2L3"
./run.sh &
cd "/data1/compbio/zlu21/AcrossTissue/ROC_Runs/Run_6LS/emb"
./run.sh &



