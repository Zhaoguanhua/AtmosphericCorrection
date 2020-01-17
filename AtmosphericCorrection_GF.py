#! usr/bin/env python
# -*- coding:utf-8 -*-
# created by zhaoguanhua 2017/09/04

import glob
import os
import sys
import tarfile            #解压缩
import json
import numpy as np
import gdal
import pdb
import math
import time
import xml.dom.minidom    #读取xml格式的影像头文件
from tqdm import tqdm     #进度条
from Py6S import *
import argparse
from base import MeanDEM

def parse_arguments(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument('--Input_dir',type=str,help='Input dir',default=None)
    parser.add_argument('--Output_dir',type=str,help='Output dir',default=None)

    return parser.parse_args(argv)

# 解压缩原始文件
def untar(fname, dirs):
    print("文件路径",fname)
    try:
        t = tarfile.open(fname)
    except Exception as e:
        print("文件%s打开失败" % fname)
    t.extractall(path=dirs)

def Block(IDataSet):
    global cols,rows,atcfiles
    #设置输出波段
    Driver = IDataSet.GetDriver()
    geoTransform1 = IDataSet.GetGeoTransform()
    ListgeoTransform1 = list(geoTransform1)
    ListgeoTransform1[5] = -ListgeoTransform1[5]
    newgeoTransform1 = tuple(ListgeoTransform1)
    proj1 = IDataSet.GetProjection()
    OutRCname = os.path.join(atcfiles,outFileName+".tif")
    outDataset = Driver.Create(OutRCname,cols,rows,4,gdal.GDT_Int32)
    outDataset.SetGeoTransform(newgeoTransform1)
    outDataset.SetProjection(proj1)
    #分别读取4个波段
    for m in range(1,5):
        ReadBand = IDataSet.GetRasterBand(m)
        outband = outDataset.GetRasterBand(m)
        outband.SetNoDataValue(-9999)
        #获取对应波段的增益gain和偏移bias
        Gain,Bias = RadiometricCalibration(m)

        #获取大气校正系数
        AtcCofa, AtcCofb, AtcCofc = AtmosphericCorrection(m)
        nBlockSize = 1024
        i = 0
        j = 0
        b = cols*rows
        #进度条参数
        XBlockcount = math.ceil(cols/nBlockSize)
        YBlockcount = math.ceil(rows/nBlockSize)
        print("第%d波段校正："%m)
        try:
            with tqdm(total=XBlockcount*YBlockcount,iterable='iterable',desc = '第%i波段:'%m) as pbar:
            # with tqdm(total=XBlockcount*YBlockcount) as pbar:
                # print(pbar)
                while i<rows:
                    while j <cols:
                        #保存分块大小
                        nXBK = nBlockSize
                        nYBK = nBlockSize

                        #最后不够分块的区域，有多少读取多少
                        if i+nBlockSize>rows:
                            nYBK = rows - i
                        if j+nBlockSize>cols:
                            nXBK=cols - j

                        #分块读取影像
                        Image = ReadBand.ReadAsArray(j,i,nXBK,nYBK)

                        outImage =np.where(Image>0,Image*Gain + Bias,-9999)

                        y = np.where(outImage!=-9999,AtcCofa * outImage - AtcCofb,-9999)
                        atcImage = np.where(y!=-9999,(y / (1 + y * AtcCofc))*10000,-9999)

                        outband.WriteArray(atcImage,j,i)

                        j=j+nXBK
                        time.sleep(1)
                        pbar.update(1)
                    j=0
                    i=i+nYBK
        except KeyboardInterrupt:
            pbar.close()
            raise
        pbar.close()

def RadiometricCalibration(BandId):
    global cols,rows,SatelliteID,SensorID,Year,ImageType,config
    if SensorID[0:3] == "WFV":
        Gain_ =config["Parameter"][SatelliteID][SensorID][Year]["gain"][BandId-1]
        Bias_ =config["Parameter"][SatelliteID][SensorID][Year]["offset"][BandId-1]
    else:
        Gain_ =config["Parameter"][SatelliteID][SensorID][Year][ImageType]["gain"][BandId-1]
        Bias_ =config["Parameter"][SatelliteID][SensorID][Year][ImageType]["offset"][BandId-1]

    return Gain_,Bias_

# 6s大气校正
def AtmosphericCorrection(BandId):
    global metedata,config,SatelliteID,SensorID
    #读取头文件
    dom = xml.dom.minidom.parse(metedata)

    # 6S模型
    s = SixS()

    # 传感器类型 自定义
    s.geometry = Geometry.User()
    s.geometry.solar_z = 90-float(dom.getElementsByTagName('SolarZenith')[0].firstChild.data)
    s.geometry.solar_a = float(dom.getElementsByTagName('SolarAzimuth')[0].firstChild.data)
    # s.geometry.view_z = float(dom.getElementsByTagName('SatelliteZenith')[0].firstChild.data)
    # s.geometry.view_a = float(dom.getElementsByTagName('SatelliteAzimuth')[0].firstChild.data)
    s.geometry.view_z = 0
    s.geometry.view_a = 0
    # 日期
    DateTimeparm = dom.getElementsByTagName('CenterTime')[0].firstChild.data
    DateTime = DateTimeparm.split(' ')
    Date = DateTime[0].split('-')
    s.geometry.month = int(Date[1])
    s.geometry.day = int(Date[2])

    # print(s.geometry)
    # 中心经纬度
    TopLeftLat = float(dom.getElementsByTagName('TopLeftLatitude')[0].firstChild.data)
    TopLeftLon = float(dom.getElementsByTagName('TopLeftLongitude')[0].firstChild.data)
    TopRightLat = float(dom.getElementsByTagName('TopRightLatitude')[0].firstChild.data)
    TopRightLon = float(dom.getElementsByTagName('TopRightLongitude')[0].firstChild.data)
    BottomRightLat = float(dom.getElementsByTagName('BottomRightLatitude')[0].firstChild.data)
    BottomRightLon = float(dom.getElementsByTagName('BottomRightLongitude')[0].firstChild.data)
    BottomLeftLat = float(dom.getElementsByTagName('BottomLeftLatitude')[0].firstChild.data)
    BottomLeftLon = float(dom.getElementsByTagName('BottomLeftLongitude')[0].firstChild.data)

    ImageCenterLat = (TopLeftLat + TopRightLat + BottomRightLat + BottomLeftLat) / 4

    # 大气模式类型
    if ImageCenterLat > -15 and ImageCenterLat < 15:
        s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.Tropical)

    if ImageCenterLat > 15 and ImageCenterLat < 45:
        if s.geometry.month > 4 and s.geometry.month < 9:
            s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeSummer)
        else:
            s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeWinter)

    if ImageCenterLat > 45 and ImageCenterLat < 60:
        if s.geometry.month > 4 and s.geometry.month < 9:
            s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.SubarcticSummer)
        else:
            s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.SubarcticWinter)

    # 气溶胶类型大陆
    s.aero_profile = AtmosProfile.PredefinedType(AeroProfile.Continental)

    # 下垫面类型
    s.ground_reflectance = GroundReflectance.HomogeneousLambertian(0.36)

    # 550nm气溶胶光学厚度,对应能见度为40km
    s.aot550 = 0.14497

    # 通过研究去区的范围去求DEM高度。
    pointUL = dict()
    pointDR = dict()
    pointUL["lat"] = max(TopLeftLat,TopRightLat,BottomRightLat,BottomLeftLat)
    pointUL["lon"] = min(TopLeftLon,TopRightLon,BottomRightLon,BottomLeftLon)
    pointDR["lat"] = min(TopLeftLat,TopRightLat,BottomRightLat,BottomLeftLat)
    pointDR["lon"] = max(TopLeftLon,TopRightLon,BottomRightLon,BottomLeftLon)
    meanDEM = (MeanDEM(pointUL, pointDR)) * 0.001

    # 研究区海拔、卫星传感器轨道高度
    s.altitudes = Altitudes()
    s.altitudes.set_target_custom_altitude(meanDEM)
    s.altitudes.set_sensor_satellite_level()

    # 校正波段（根据波段名称）
    if BandId == 1:
        SRFband = config["Parameter"][SatelliteID][SensorID]["SRF"]["1"]
        s.wavelength = Wavelength(0.450,0.520,SRFband)

    elif BandId == 2:
        SRFband = config["Parameter"][SatelliteID][SensorID]["SRF"]["2"]

        s.wavelength = Wavelength(0.520,0.590,SRFband)

    elif BandId == 3:
        SRFband = config["Parameter"][SatelliteID][SensorID]["SRF"]["3"]

        s.wavelength = Wavelength(0.630,0.690,SRFband)

    elif BandId == 4:
        SRFband = config["Parameter"][SatelliteID][SensorID]["SRF"]["4"]
        s.wavelength = Wavelength(0.770,0.890,SRFband)

    s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromReflectance(-0.1)

    # 运行6s大气模型
    s.run()
    xa = s.outputs.coef_xa
    xb = s.outputs.coef_xb
    xc = s.outputs.coef_xc
    # x = s.outputs.values
    return (xa, xb, xc)

if __name__ == '__main__':

    #读取辐射校正和大气校正所需参数:增益、偏移和光谱响应函数
    config = json.load(open("RadiometricCorrectionParameter.json"))

    #输入数据路径
    InputFilePath = parse_arguments(sys.argv[1:]).Input_dir
    OutputFilePath = parse_arguments(sys.argv[2:]).Output_dir

    for root,dirs,tarFiles in os.walk(InputFilePath):
        pass

    for tarFile in tarFiles:
        print(tarFile)
        filename = os.path.basename(tarFile)
        fileType = filename[0:2]
        if fileType == 'GF':
            GFType = filename[4:7]
            intputname = os.path.join(InputFilePath,filename)
            outFileName = filename[:-7]
            outname = os.path.join(InputFilePath,outFileName)
            atcfiles = os.path.join(OutputFilePath,outFileName)

            print("文件"+filename+"开始解压缩")
            if GFType == 'WFV':
                try:
                    untar(intputname,outname)
                except Exception as e:
                    continue
                tiffFile = glob.glob(outname + "/*.tiff")[0]
                metedata = glob.glob(outname+"/*.xml")[0]

            elif GFType == 'PMS':
                try:
                    untar(intputname, outname)
                except Exception as e:
                    pass

                tiffFile = glob.glob(outname + "/*MSS*.tiff")[0]
                metedata = glob.glob(outname+"/*MSS*.xml")[0]

            try:
                os.mkdir(atcfiles)
            except Exception as e:
                pass
            print(filename+"解压缩完成")

            try:
                IDataSet = gdal.Open(tiffFile)
            except Exception as e:
                print("文件%S打开失败" % tiffFile)

            cols = IDataSet.RasterXSize
            rows = IDataSet.RasterYSize

            SatelliteID = filename[0:3]
            SensorID = filename[4:8]
            Year = filename[22:26]
            ImageType =os.path.basename(tiffFile)[-9:-6]

            Block(IDataSet)

