# RNAseqFlow

## 描述 

## 安装

```R
devtools::install_github("Feng-Zhang/RNAseqFlow") #安装稳定版
devtools::install_github("Feng-Zhang/RNAseqFlow@dev") #安装开发版
```



## 使用



## to do

- RNAseqFlow要安装的软件太多了，需要减少依赖包数量
- 添加batch相关信息，将RNAseqWorkflowDEseq2的内容合并进来。分支：feature_20210702_batch
- 修复Warning messages。 当运行DEseqObj时会出现:  
```
1: In rownames(col_data) == colnames(count_data) :
  长的对象长度不是短的对象长度的整倍数
```
