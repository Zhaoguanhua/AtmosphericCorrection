#! usr/bin/env python
# -*- coding:utf-8 -*-
# created by zhaoguanhua 2017/9/25
# AtmosphericCorrection for Landsat8

import glob
import os
import sys
import tarfile
import re
import gdal
import numpy
from Py6S import *
from osgeo import gdal
import pdb
import shutil

# 解压缩原始文件
def untar(fname, dirs):
    try:
        t = tarfile.open(fname)
    except Exception as e:
        print("文件%s打开失败" % fname)
    t.extractall(path=dirs)

# 逐波段辐射定标
def RadiometricCalibration(BandId):
    # LandSat8 TM辐射定标参数
    global data2,ImgRasterData
    parameter_OLI = numpy.zeros((11,2))

    #计算辐射亮度参数
    # parameter_OLI[0,0] = float(''.join(re.findall('RADIANCE_MULT_BAND_1.+',data2)).split("=")[1])
    parameter_OLI[1,0] = float(''.join(re.findall('RADIANCE_MULT_BAND_2.+',data2)).split("=")[1])
    parameter_OLI[2,0] = float(''.join(re.findall('RADIANCE_MULT_BAND_3.+',data2)).split("=")[1])
    parameter_OLI[3,0] = float(''.join(re.findall('RADIANCE_MULT_BAND_4.+',data2)).split("=")[1])
    parameter_OLI[4,0] = float(''.join(re.findall('RADIANCE_MULT_BAND_5.+',data2)).split("=")[1])
    parameter_OLI[5,0] = float(''.join(re.findall('RADIANCE_MULT_BAND_6.+',data2)).split("=")[1])
    parameter_OLI[6,0] = float(''.join(re.findall('RADIANCE_MULT_BAND_7.+',data2)).split("=")[1])
    parameter_OLI[7,0] = float(''.join(re.findall('RADIANCE_MULT_BAND_8.+',data2)).split("=")[1])
    parameter_OLI[8,0] = float(''.join(re.findall('RADIANCE_MULT_BAND_9.+',data2)).split("=")[1])
    parameter_OLI[9,0] = float(''.join(re.findall('RADIANCE_MULT_BAND_10.+',data2)).split("=")[1])
    parameter_OLI[10,0] = float(''.join(re.findall('RADIANCE_MULT_BAND_11.+',data2)).split("=")[1])


    # parameter_OLI[0,1] = float(''.join(re.findall('RADIANCE_ADD_BAND_1.+',data2)).split("=")[1])
    parameter_OLI[1,1] = float(''.join(re.findall('RADIANCE_ADD_BAND_2.+',data2)).split("=")[1])
    parameter_OLI[2,1] = float(''.join(re.findall('RADIANCE_ADD_BAND_3.+',data2)).split("=")[1])
    parameter_OLI[3,1] = float(''.join(re.findall('RADIANCE_ADD_BAND_4.+',data2)).split("=")[1])
    parameter_OLI[4,1] = float(''.join(re.findall('RADIANCE_ADD_BAND_5.+',data2)).split("=")[1])
    parameter_OLI[5,1] = float(''.join(re.findall('RADIANCE_ADD_BAND_6.+',data2)).split("=")[1])
    parameter_OLI[6,1] = float(''.join(re.findall('RADIANCE_ADD_BAND_7.+',data2)).split("=")[1])
    parameter_OLI[7,1] = float(''.join(re.findall('RADIANCE_ADD_BAND_8.+',data2)).split("=")[1])
    parameter_OLI[8,1] = float(''.join(re.findall('RADIANCE_ADD_BAND_9.+',data2)).split("=")[1])
    parameter_OLI[9,1] = float(''.join(re.findall('RADIANCE_ADD_BAND_10.+',data2)).split("=")[1])
    parameter_OLI[10,1] = float(''.join(re.findall('RADIANCE_ADD_BAND_11.+',data2)).split("=")[1])

    if len(BandId) ==8:
        n = int(BandId[2])
    else:
        n = int(BandId[1:3])
    Gain = parameter_OLI[n - 1,0]
    Bias = parameter_OLI[n - 1,1]

    RaCal = numpy.where(ImgRasterData>0 ,Gain * ImgRasterData + Bias,-9999)
    return (RaCal)

#计算表观反射率
def TOAReflectance(BandId):
    # LandSat8 TM辐射定标参数
    global data2,ImgRasterData
    parameter_OLI = numpy.zeros((9,2))

    #计算表观反射率参数
    parameter_OLI[0,0] = float(''.join(re.findall('REFLECTANCE_MULT_BAND_1.+',data2)).split("=")[1])
    parameter_OLI[1,0] = float(''.join(re.findall('REFLECTANCE_MULT_BAND_2.+',data2)).split("=")[1])
    parameter_OLI[2,0] = float(''.join(re.findall('REFLECTANCE_MULT_BAND_3.+',data2)).split("=")[1])
    parameter_OLI[3,0] = float(''.join(re.findall('REFLECTANCE_MULT_BAND_4.+',data2)).split("=")[1])
    parameter_OLI[4,0] = float(''.join(re.findall('REFLECTANCE_MULT_BAND_5.+',data2)).split("=")[1])
    parameter_OLI[5,0] = float(''.join(re.findall('REFLECTANCE_MULT_BAND_6.+',data2)).split("=")[1])
    parameter_OLI[6,0] = float(''.join(re.findall('REFLECTANCE_MULT_BAND_7.+',data2)).split("=")[1])
    parameter_OLI[7,0] = float(''.join(re.findall('REFLECTANCE_MULT_BAND_8.+',data2)).split("=")[1])
    parameter_OLI[8,0] = float(''.join(re.findall('REFLECTANCE_MULT_BAND_9.+',data2)).split("=")[1])

    parameter_OLI[0,1] = float(''.join(re.findall('REFLECTANCE_ADD_BAND_1.+',data2)).split("=")[1])
    parameter_OLI[1,1] = float(''.join(re.findall('REFLECTANCE_ADD_BAND_2.+',data2)).split("=")[1])
    parameter_OLI[2,1] = float(''.join(re.findall('REFLECTANCE_ADD_BAND_3.+',data2)).split("=")[1])
    parameter_OLI[3,1] = float(''.join(re.findall('REFLECTANCE_ADD_BAND_4.+',data2)).split("=")[1])
    parameter_OLI[4,1] = float(''.join(re.findall('REFLECTANCE_ADD_BAND_5.+',data2)).split("=")[1])
    parameter_OLI[5,1] = float(''.join(re.findall('REFLECTANCE_ADD_BAND_6.+',data2)).split("=")[1])
    parameter_OLI[6,1] = float(''.join(re.findall('REFLECTANCE_ADD_BAND_7.+',data2)).split("=")[1])
    parameter_OLI[7,1] = float(''.join(re.findall('REFLECTANCE_ADD_BAND_8.+',data2)).split("=")[1])
    parameter_OLI[8,1] = float(''.join(re.findall('REFLECTANCE_ADD_BAND_9.+',data2)).split("=")[1])

    n = int(BandId[1])
    Gain = parameter_OLI[n - 1,0]
    Bias = parameter_OLI[n - 1,1]


    SunElevationFactor = numpy.sin(float(''.join(re.findall('SUN_ELEVATION.+',data2)).split("=")[1])/180*numpy.pi)

    TOARef = numpy.where(ImgRasterData>0,(Gain * ImgRasterData + Bias)/SunElevationFactor,-9999)

    return (TOARef)

# 6s大气校正
def AtmosphericCorrection(BandId):
    global data
    # 6S模型
    s = SixS()

    s.geometry = Geometry.User()
    s.geometry.solar_z = 90-float(''.join(re.findall('SUN_ELEVATION.+',data2)).split("=")[1])
    s.geometry.solar_a = float(''.join(re.findall('SUN_AZIMUTH.+',data2)).split("=")[1])
    s.geometry.view_z = 0
    s.geometry.view_a = 0


    # 日期
    Dateparm = ''.join(re.findall('DATE_ACQUIRED.+',data2)).split("=")
    Date = Dateparm[1].split('-')

    s.geometry.month = int(Date[1])
    s.geometry.day = int(Date[2])

    # 中心经纬度
    point1lat = float(''.join(re.findall('CORNER_UL_LAT_PRODUCT.+',data2)).split("=")[1])
    point1lon = float(''.join(re.findall('CORNER_UL_LON_PRODUCT.+',data2)).split("=")[1])
    point2lat = float(''.join(re.findall('CORNER_UR_LAT_PRODUCT.+',data2)).split("=")[1])
    point2lon = float(''.join(re.findall('CORNER_UR_LON_PRODUCT.+',data2)).split("=")[1])
    point3lat = float(''.join(re.findall('CORNER_LL_LAT_PRODUCT.+',data2)).split("=")[1])
    point3lon = float(''.join(re.findall('CORNER_LL_LON_PRODUCT.+',data2)).split("=")[1])
    point4lat = float(''.join(re.findall('CORNER_LR_LAT_PRODUCT.+',data2)).split("=")[1])
    point4lon = float(''.join(re.findall('CORNER_LR_LON_PRODUCT.+',data2)).split("=")[1])

    sLongitude = (point1lon + point2lon + point3lon + point4lon) / 4
    sLatitude = (point1lat + point2lat + point3lat + point4lat) / 4

    # 大气模式类型
    if sLatitude > -15 and sLatitude <= 15:
        s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.Tropical)

    if sLatitude > 15 and sLatitude <= 45:
        if s.geometry.month > 4 and s.geometry.month <= 9:
            s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeSummer)
        else:
            s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeWinter)

    if sLatitude > 45 and sLatitude <= 60:
        if s.geometry.month > 4 and s.geometry.month <= 9:
            s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.SubarcticSummer)
        else:
            s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.SubarcticWinter)

    # 气溶胶类型大陆
    s.aero_profile = AtmosProfile.PredefinedType(AeroProfile.Continental)

    # 目标地物？？？？？？
    s.ground_reflectance = GroundReflectance.HomogeneousLambertian(0.36)

    # 550nm气溶胶光学厚度,根据日期从MODIS处获取。
    #s.visibility=40.0
    s.aot550 = 0.14497

    # 通过研究去区的范围去求DEM高度。
    pointUL = dict()
    pointDR = dict()
    pointUL["lat"] = point1lat
    pointUL["lon"] = point1lon
    pointDR["lat"] = point4lat
    pointDR["lon"] = point2lon
    meanDEM = (MeanDEM(pointUL, pointDR)) * 0.001

    # 研究区海拔、卫星传感器轨道高度
    s.altitudes = Altitudes()
    s.altitudes.set_target_custom_altitude(meanDEM)
    s.altitudes.set_sensor_satellite_level()

    # 校正波段（根据波段名称）
    if BandId == 'B1.TIF':
        s.wavelength = Wavelength(PredefinedWavelengths.LANDSAT_OLI_B1)

    elif BandId == 'B2.TIF':
        s.wavelength = Wavelength(PredefinedWavelengths.LANDSAT_OLI_B2)

    elif BandId == 'B03.tiff':
        s.wavelength = Wavelength(PredefinedWavelengths.LANDSAT_OLI_B3)

    elif BandId == 'B04.tiff':
        s.wavelength = Wavelength(PredefinedWavelengths.LANDSAT_OLI_B4)

    elif BandId == 'B05.tiff':
        s.wavelength = Wavelength(PredefinedWavelengths.LANDSAT_OLI_B5)

    elif BandId == 'B6.TIF':
        s.wavelength = Wavelength(PredefinedWavelengths.LANDSAT_OLI_B6)

    elif BandId == 'B7.TIF':
        s.wavelength = Wavelength(PredefinedWavelengths.LANDSAT_OLI_B7)

    elif BandId == 'B8.TIF':
        s.wavelength = Wavelength(PredefinedWavelengths.LANDSAT_OLI_B8)

    elif BandId == 'B9.TIF':
        s.wavelength = Wavelength(PredefinedWavelengths.LANDSAT_OLI_B9)

    # 下垫面非均一、朗伯体
    s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromReflectance(-0.1)

    # 运行6s大气模型
    s.run()

    xa = s.outputs.coef_xa
    xb = s.outputs.coef_xb
    xc = s.outputs.coef_xc
    x = s.outputs.values
    print(x)
    return (xa, xb, xc)

def MeanDEM(pointUL, pointDR):
    # 打开DEM数据
    try:
        DEMIDataSet = gdal.Open("GMTED2km.tif")
    except Exception as e:
        pass

    DEMBand = DEMIDataSet.GetRasterBand(1)
    cols = DEMIDataSet.RasterXSize
    rows = DEMIDataSet.RasterYSize

    geotransform = DEMIDataSet.GetGeoTransform()
    # DEM分辨率
    pixelWidth = geotransform[1]
    pixelHight = geotransform[5]

    # DEM起始点：左上角，X：经度，Y：纬度
    originX = geotransform[0]
    originY = geotransform[3]

    # 研究区左上角在矩阵中的位置
    yoffset1 = int((originY - pointUL['lat']) / pixelWidth)
    xoffset1 = int((pointUL['lon'] - originX) / (-pixelHight))

    # 研究区右下角在矩阵中的位置
    yoffset2 = int((originY - pointDR['lat']) / pixelWidth)
    xoffset2 = int((pointDR['lon'] - originX) / (-pixelHight))

    # 研究区矩阵行列数
    xx = xoffset2 - xoffset1
    yy = yoffset2 - yoffset1


    # 读取研究区内的数据，并计算高程
    DEMRasterData = DEMBand.ReadAsArray(xoffset1, yoffset1, xx, yy)

    MeanAltitude = numpy.mean(DEMRasterData)
    return MeanAltitude

def CloudMaskScore():

    mask = BrightTemp == -9999

    #得分1：Blue
    BluePart = numpy.ma.array((TOARefRasterBlue-0.1)/0.2,mask=mask)
    BluePart.fill_value=-9999
    ScorePart1 = numpy.where(BluePart.filled()>1,1,BluePart.filled())

    #得分2：Red、Blue、Green
    RGBPart = numpy.ma.array((TOARefRasterBlue+TOARefRasterGreen+TOARefRasterRed-0.2)/0.6,mask=mask)
    RGBPart.fill_value=-9999
    ScorePart2 = numpy.where(RGBPart.filled()>ScorePart1,ScorePart1,RGBPart.filled())

    #得分3：Nir、Swir1、Swir2
    NSSPart = numpy.ma.array((TOARefRasterNir+TOARefRasterSwir1+TOARefRasterSwir2-0.3)/0.5,mask=mask)
    NSSPart.fill_value=-9999
    ScorePart3 = numpy.where(NSSPart.filled()>ScorePart2,ScorePart2,NSSPart.filled())

    #得分4:temperature
    TempPart = numpy.ma.array((BrightTemp-300)/(-10),mask=mask)
    TempPart.fill_value=-9999
    ScorePart4 = numpy.where(TempPart.filled()>ScorePart3,ScorePart3,TempPart.filled())

    #得分5NDSI：Green、TOARefRasterSwir1
    NDSIPart1 = numpy.ma.array((TOARefRasterGreen- TOARefRasterSwir1)/(TOARefRasterGreen+TOARefRasterSwir1),mask=mask)
    NDSIPart2 = numpy.ma.array((NDSIPart1-0.8)/(-0.2),mask=mask)
    NDSIPart2.fill_value=-9999
    ScorePart5 = numpy.where(NDSIPart2.filled()>ScorePart4,ScorePart4,NDSIPart2.filled())

    ScoreCloud = numpy.where(ScorePart5!=-9999,1- ScorePart5,-9999)
    return ScoreCloud

if __name__ == '__main__':

    #输入数据路径
    RootOutName = sys.argv[2]
    RootInputPath = sys.argv[1]

    Contronal=0
    #创建日志文件
    LogFile = open(os.path.join(RootOutName,'log.text'),'w')

    for root,dirs,RSFiles in os.walk(RootInputPath):

        #判断是否进入最底层
        if len(dirs)==0:
            #根据输入输出路径建立生成新文件的路径
            RootInputPathList = RootInputPath.split('/')
            RootList = root.split('/')
            StartList = len(RootInputPathList)
            EndList = len(RootList)
            outname = RootOutName
            for i in range(StartList,EndList):
                if os.path.exists(os.path.join(outname,RootList[i]))==False:
                    os.makedirs(os.path.join(outname,RootList[i]))
                    outname=os.path.join(outname,RootList[i])
                else:
                    outname=os.path.join(outname,RootList[i])

            #判断文件是否都存在
            CloudScoreFile = os.path.join(outname,RootList[-1]+'_CloudScore.TIF')

            if os.path.isfile(CloudScoreFile):
                print(root+'计算完成')
                continue
            else:
                MeteData = os.path.join(root,'MTL.txt')
                f = open(MeteData)
                data = f.readlines()
                data2 =' '.join(data)
                
                shutil.copyfile(MeteData,os.path.join(outname,RootList[-1]+'MTL.txt'))

                for tifFile in RSFiles:
                    # print(tifFile)
                    if tifFile[-5:] == '.tiff':
                        BandId = (os.path.basename(tifFile))
                        # print(BandId)

                        #捕捉打开数据出错异常
                        try:
                            IDataSet = gdal.Open(os.path.join(root,tifFile))
                        except Exception as e:
                            print("文件%S打开失败" % tifFile)
                            LogFile.write('\n'+os.path.join(root,tifFile)+'数据打开失败')
                        
                        if IDataSet == None:
                            LogFile.write('\n'+os.path.join(root,tifFile)+'数据集读取为空')
                            continue
                        else:
                            #获取行列号
                            cols = IDataSet.RasterXSize
                            rows = IDataSet.RasterYSize
                            ImgBand = IDataSet.GetRasterBand(1)
                            ImgRasterData = ImgBand.ReadAsArray(0, 0, cols, rows)

                            if ImgRasterData is None:
                                LogFile.write('\n'+os.path.join(root,tifFile)+'栅格数据为空')
                                continue
                            else:
                                if BandId =='B02.tiff':
                                    # TOARefRasterBlue = TOAReflectance(BandId)
                                    RaCalRaster = RadiometricCalibration(BandId)
                                    Contronal = Contronal + 1
                                elif BandId =='B03.tiff':
                                    # TOARefRasterGreen = TOAReflectance(BandId)
                                    RaCalRaster = RadiometricCalibration(BandId)
                                    Contronal = Contronal + 1
                                elif BandId =='B04.tiff':
                                    # TOARefRasterRed = TOAReflectance(BandId)
                                    RaCalRaster = RadiometricCalibration(BandId)
                                    Contronal = Contronal + 1
                                elif BandId =='B05.tiff':
                                    # TOARefRasterNir = TOAReflectance(BandId)
                                    RaCalRaster = RadiometricCalibration(BandId)
                                    Contronal = Contronal + 1
                                # elif BandId =='B6.TIF':
                                #     TOARefRasterSwir1 = TOAReflectance(BandId)
                                #     Contronal = Contronal + 1
                                # elif BandId =='B7.TIF':
                                #     TOARefRasterSwir2 = TOAReflectance(BandId)
                                #     Contronal = Contronal + 1
                                # elif tifFile[-7:] =='B10.TIF':
                                #     RaCalRaster = RadiometricCalibration(BandId)
                                #     Contronal = Contronal + 1
                                #     BrightTemp = numpy.where(RaCalRaster!=-9999,1321.08/numpy.log(774.89/RaCalRaster+1),-9999)
                                    # print("亮温计算完成")

                                if BandId == 'B02.tiff'or BandId == 'B03.tiff'or BandId == 'B04.tiff'or BandId == 'B05.tiff':
                                    #设置输出文件路径
                                    outFilename=os.path.join(outname,os.path.basename(tifFile))

                                    #如果文件存在就跳过，进行下一波段操作
                                    if os.path.isfile(outFilename):
                                        print("%s已经完成" % outFilename)
                                        continue
                                    else:
                                        # #辐射校正
                                        # RaCalRaster = RadiometricCalibration(tifFile, BandId)
                                        #大气校正
                                        a, b, c = AtmosphericCorrection(BandId)
                                        y = numpy.where(RaCalRaster!=-9999,a * RaCalRaster - b,-9999)
                                        atc = numpy.where(y!=-9999,(y / (1 + y * c))*10000,-9999)
                                        
                                        driver = IDataSet.GetDriver()
                                        #输出栅格数据集
                                        outDataset = driver.Create(outFilename, cols, rows, 1, gdal.GDT_Int16)

                                        # 设置投影信息，与原数据一样
                                        geoTransform = IDataSet.GetGeoTransform()
                                        outDataset.SetGeoTransform(geoTransform)
                                        proj = IDataSet.GetProjection()
                                        outDataset.SetProjection(proj)

                                        outband = outDataset.GetRasterBand(1)
                                        outband.SetNoDataValue(-9999)
                                        outband.WriteArray(atc, 0, 0)
                print(root+'计算完成')

                # if Contronal == 7:
                #     print(Contronal)
                #     #设置输出文件路径
                #     outFilename=os.path.join(outname,os.path.basename(tifFile)[0:41]+'CloudScore.TIF')
                #     CloudScoreFile = os.path.join(outname,RootList[-1]+'_CloudScore.TIF')

                #     CloudScoreData = CloudMaskScore()
                #     driver = IDataSet.GetDriver()
                #     #输出栅格数据集
                #     CloudDataset = driver.Create(CloudScoreFile, cols, rows, 1, gdal.GDT_Float32)

                #     # 设置投影信息，与原数据一样
                #     geoTransform = IDataSet.GetGeoTransform()
                #     proj = IDataSet.GetProjection()
                #     CloudDataset.SetGeoTransform(geoTransform)
                #     CloudDataset.SetProjection(proj)

                #     outband = CloudDataset.GetRasterBand(1)
                #     outband.SetNoDataValue(-9999)
                #     outband.WriteArray(CloudScoreData, 0, 0)
                #     print('影像'+outFilename + '处理完成')
                #     RasterData = None
                #     Contronal=0
                # else:
                #     Contronal=0

    #关闭日志文件
    LogFile.close()


