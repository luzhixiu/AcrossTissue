
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



cd "/data1/compbio/zlu21/AcrossTissue/ROC_Runs/Run_L2L3_Combined/01_24_No_Lower_Bound/adult"
./run.sh &
cd "/data1/compbio/zlu21/AcrossTissue/ROC_Runs/Run_L2L3_Combined/01_24_No_Lower_Bound/dauer"
./run.sh &
cd  "/data1/compbio/zlu21/AcrossTissue/ROC_Runs/Run_L2L3_Combined/01_24_No_Lower_Bound/L1"
./run.sh &
cd  "/data1/compbio/zlu21/AcrossTissue/ROC_Runs/Run_L2L3_Combined/01_24_No_Lower_Bound/L2L3"
./run.sh &
cd "/data1/compbio/zlu21/AcrossTissue/ROC_Runs/Run_L2L3_Combined/01_24_No_Lower_Bound/emb"
./run.sh &



