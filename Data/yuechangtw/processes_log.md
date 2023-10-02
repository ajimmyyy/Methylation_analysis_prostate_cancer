# 完整流程
1. 開發環境建置
2. 在RStudio輸入IDAT檔（分為測試集test和訓練集train），經ChAMP套件之norm正規化後輸出存有甲基化後beta值的csv檔。以及用ChAMP.DMP分析輸出各位點資料之csv檔。
myLoad <- champ.load("PATH/train")#, SampleCutoff = 0.2)#, arraytype = "EPIC")
myNorm <- champ.norm(beta=myLoad$beta,plotBMIQ=FALSE,cores=5)#, arraytype = "EPIC")
write.csv(myNorm,file="all_beta_normalized.csv",quote=F,row.names=T)
myDMP <- champ.DMP(beta=myNorm, pheno=myLoad$pd$Sample_Group)#, arraytype = "EPIC")
write.csv(myDMP[[1]], file="DMP_result_TN.csv", quote=F)

3. 
