# 基于6s模型的影像大气校正工程
## OVERVIEW
调用py6s接口，自动读取影像头文件信息，对遥感影像进行大气校正批处理。
## 环境 & 依赖
conda create -n py6s-env -c conda-forge py6s

## 脚本说明

* AtmosphericCorrection.py 针对landsat8影像,已经可以工程化使用。
* AtmosphericCorrection_Sentinel.py 针对Sentinel影像，已经可以工程化使用。
* AtmosphericCorrection_GF.py 针对GF1、2影像，已经可以工程化使用。

## 测试

```
source activate py6s-env
python .../AtmosphericCorrection/AtmosphericCorrection_Lansat8.py Input_dir=输入路径 Output_dir=输出路径
python .../AtmosphericCorrection/AtmosphericCorrection_Sentinel.py Input_dir=输入路径 Output_dir=输出路径
python .../AtmosphericCorrection/AtmosphericCorrection_GF.py Input_dir=输入路径 Output_dir=输出路径
```

