## Usage

![](usage.png)

python ./hed.py -t hla.result.tsv –s {sample} –f ./db/db.tsv –d ./bin/AAdistMatrix_Grantham.cnv -o ./Output

## result

```
Output
├── hla_type_change.txt
├── {sample}.HLA.fa
├── {sample}.HLA.fa_PairwiseDistanceMatrix.txt
├── {sample}.HLA.fa_PairwiseDistanceList.txt
├── {sample}.HLA.fa_PairwiseDistanceList.txt.mod
└── {sample}.hed.tsv
```

### cat

`{sample}.hed.tsv`

HLA | A1      | A2      | B1      | B2      | C1      | C2      | Reads | Objective | QC            | HED-A | HED-B | HED-C | HED-TOTAL
--- | ------- | ------- | ------- | ------- | ------- | ------- | ----- | --------- | ------------- | ----- | ----- | ----- | ---------
0   | A*02:06 | A*24:20 | B*46:01 | B*15:11 | C*03:03 | C*01:02 | 18745 | 18070.18  | Qualification | 10.73 | 3.31  | 5.92  | 6.65

