#!usr/bin/env python
# -*- coding: utf-8 -*-
"""
@Author  : zhaoguanhua
@Email   :
@Time    : 2021/4/2 16:01
@File    : AtmosphericCorrection_multiprocess.py
@Software: PyCharm
"""
from multiprocessing import Pool
import os,time,random
import glob
import json
import gdal
from AtmosphericCorrection_GF import untar,Block

def long_time_task(name):
    print('run task {} {} ...'.format(name,os.getpid()))
    start=time.time()
    time.sleep(random.random()*3)
    end = time.time()
    print('Task {} runs {} seconds'.format(name,(end-start)))

def atmospheric_correction(file_path,input_dir,output_dir,config):
    print('run task {} {} ...'.format(file_path, os.getpid()))

    file_name = os.path.basename(file_path)
    fileType = file_name[0:2]
    filename_split = file_name.split("_")

    GFType = filename_split[1][:3]     #传感器
    untar_dirname = file_name[:-7]                 #解压后影像文件夹名
    untar_dir = os.path.join(input_dir, untar_dirname)#解压后影像文件夹路径

    print("文件" + file_path + "开始解压缩")

    try:
        untar(file_path, untar_dir)
    except Exception as e:
        pass

    if GFType == 'WFV':
        tiffFile = glob.glob(os.path.join(untar_dir, "*.tiff"))[0]
        metedata = glob.glob(os.path.join(untar_dir, "*.xml"))[0]

    elif GFType == 'PMS':
        tiffFile = glob.glob(os.path.join(untar_dir, "*MSS*.tiff"))[0]
        metedata = glob.glob(os.path.join(untar_dir, "*MSS*.xml"))[0]

    atcfile_dir = os.path.join(output_dir,untar_dirname) #大气校正结果文件夹

    try:
        os.mkdir(atcfile_dir)
    except Exception as e:
        pass
    # print(filename + "解压缩完成")

    try:
        IDataSet = gdal.Open(tiffFile)
        print(IDataSet)
    except Exception as e:
        print("文件%S打开失败" % tiffFile)

    # Block(IDataSet)
    ImageType = os.path.basename(tiffFile)[-9:-6]

    Block(IDataSet, filename_split, atcfile_dir, ImageType, config, metedata,untar_dirname)


if __name__ == '__main__':
    print('Parent process {}'.format(os.getpid()))

    script_path = os.path.split(os.path.realpath(__file__))[0]
    #读取辐射校正和大气校正所需参数:增益、偏移和光谱响应函数
    config_file = os.path.join(script_path,"RadiometricCorrectionParameter.json")
    config = json.load(open(config_file))

    input_dir=r"D:\test_data"
    output_dir=r"D:\temp"

    #获取影像列表
    GF_files= glob.glob(os.path.join(input_dir,"*.tar.gz"))
    print(GF_files)

    #进程池
    p=Pool(2)
    for gf_file_path in GF_files:
        p.apply_async(atmospheric_correction,args=(gf_file_path,input_dir,output_dir,config,))


    print('Waiting for all subprocesses done...')
    p.close()
    p.join()
    print("All subprocesses done")