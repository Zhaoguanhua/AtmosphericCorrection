# 基于6s模型的影像大气校正工程
## OVERVIEW
调用py6s接口，自动读取影像头文件信息，对遥感影像进行大气校正批处理。
## 环境 & 依赖
python版本3.6  
conda install gdal  
conda install -c conda-forge py6s


## 脚本说明

* AtmosphericCorrection_Landsat8.py 针对landsat8影像,已经可以工程化使用。
* AtmosphericCorrection_Sentinel.py 针对Sentinel影像，已经可以工程化使用。
* AtmosphericCorrection_GF.py 针对GF1、2影像，已经可以工程化使用。   
为了减少校正结果存储空间，程序中将大气校正的结果放大了10000倍。


## 测试

```
python .../AtmosphericCorrection/AtmosphericCorrection_Lansat8.py Input_dir=输入路径 Output_dir=输出路径
python .../AtmosphericCorrection/AtmosphericCorrection_Sentinel.py Input_dir=输入路径 Output_dir=输出路径
python .../AtmosphericCorrection/AtmosphericCorrection_GF.py Input_dir=输入路径 Output_dir=输出路径
```
更新说明：针对标准存储格式的sentinel-2大气校正做了修改  
1、文件存储结构图  
<img src="https://github.com/Zhaoguanhua/AtmosphericCorrection/blob/master/img/sentinel-2_ImageTree.png" width=50%>  

2、参数使用  
输入参数中的输入路径需要到L1C_T51TUE_A004877_20180211T025320这一级(与压缩文件同名)，这一级下的IMG_DATA文件夹存储了
各波段影像文件，MTD_TL.xml是影像元文件  
<img src="https://github.com/Zhaoguanhua/AtmosphericCorrection/blob/master/img/sentinel-2_AC.png" width=50%>
## <font color=red>注意</font>
直接在pycharm测试可能会有bug，建议windows用户直接在conda自带的Anaconda Prompt工具中测试，mac可直接在终端里测试。
