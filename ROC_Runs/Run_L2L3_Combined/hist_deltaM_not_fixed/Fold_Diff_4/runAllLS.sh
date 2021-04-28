
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


curDir="/data1/compbio/zlu21/AcrossTissue/ROC_Runs/Run_L2L3_Combined/Fold_Diff_4"

cd "${curDir}/adult"
./run.sh &
cd "${curDir}/dauer"
./run.sh &
cd  "${curDir}/L1"
./run.sh &
cd  "${curDir}/L2L3"
./run.sh &
cd "${curDir}/emb"
./run.sh &



