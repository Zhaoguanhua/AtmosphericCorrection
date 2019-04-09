#! usr/bin/env python
# -*- coding:utf-8 -*-
# created by zhaoguanhua 2017/10/9
# AtmosphericCorrection for Sentinel-2A

import glob
import os
import sys
import tarfile
import re
import numpy
from Py6S import *
from osgeo import gdal
from osgeo import osr
import xml.dom.minidom
import pdb
import shutil

def TOAReflectanceToTOARadiance(BandId):
    '''
    将表观反射率转换为表观辐射亮度
    '''
    global dom
    #太阳辐照度
    EsBand = numpy.zeros((14))
    EsBand[1] = 1913.57
    EsBand[2] = 1941.63
    EsBand[3] = 1822.61
    EsBand[4] = 1512.79
    EsBand[5] = 1425.56
    EsBand[6] = 1288.32
    EsBand[7] = 1163.19
    EsBand[8] = 1036.39
    EsBand[9] = 955.19
    EsBand[10]= 813.04 
    EsBand[11]= 367.15
    EsBand[12]= 245.59
    EsBand[13]= 85.25

    #日地相对距离，在1左右波动，暂时用1代替
    #公式：D=1+0.0167*sin(2*pi*(days-39.5)/360)  days是儒略日
    D = 1
    #太阳天顶角
    SunZenith = float(dom.getElementsByTagName('ZENITH_ANGLE')[0].firstChild.data)
    Us = numpy.cos(SunZenith/180*numpy.pi)

    Radiance =numpy.where(ImgRasterData!=0,(ImgRasterData/10000)*Us*EsBand[BandId]/(D*D*numpy.pi),-9999)
    return Radiance

# 6s大气校正
def AtmosphericCorrection(BandId):
    '''
    调用6s模型，给各参数赋值，得到大气校正参数
    '''
    global dom,SixsInputParameter
    # 6S模型
    s = SixS()

    s.geometry = Geometry.User()
    s.geometry.solar_z = SixsInputParameter["SolarZenithAngle"]
    s.geometry.solar_a = SixsInputParameter["SolarAzimuthAngle"]
    s.geometry.view_z = SixsInputParameter["ViewZenithAngle"][BandId]
    s.geometry.view_a = SixsInputParameter["ViewAzimuthAngle"][BandId]

    # 日期:月、日
    s.geometry.month = SixsInputParameter["ImgMonth"]
    s.geometry.day = SixsInputParameter["ImgDay"]

    #大气模式类型
    s.atmos_profile = AtmosProfile.PredefinedType(SixsInputParameter["AtmosphericProfile"])

    # 气溶胶类型大陆
    s.aero_profile = AtmosProfile.PredefinedType(AeroProfile.Continental)

    # 目标地物
    s.ground_reflectance = GroundReflectance.HomogeneousLambertian(0.36)

    # 550nm气溶胶光学厚度,根据日期从MODIS处获取。
    #s.visibility=40.0
    s.aot550 = 0.14497

    # 研究区海拔、卫星传感器轨道高度
    s.altitudes = Altitudes()
    s.altitudes.set_target_custom_altitude(SixsInputParameter["meanDEM"])
    s.altitudes.set_sensor_satellite_level()

    # 校正波段（根据波段名称）
    if BandId == 1:
        s.wavelength = Wavelength(PredefinedWavelengths.S2A_MSI_01)

    elif BandId == 2:
        s.wavelength = Wavelength(PredefinedWavelengths.S2A_MSI_02)

    elif BandId == 3:
        s.wavelength = Wavelength(PredefinedWavelengths.S2A_MSI_03)

    elif BandId == 4:
        s.wavelength = Wavelength(PredefinedWavelengths.S2A_MSI_04)

    elif BandId == 5:
        s.wavelength = Wavelength(PredefinedWavelengths.S2A_MSI_05)

    elif BandId == 6:
        s.wavelength = Wavelength(PredefinedWavelengths.S2A_MSI_06)

    elif BandId == 7:
        s.wavelength = Wavelength(PredefinedWavelengths.S2A_MSI_07)

    elif BandId == 8:
        s.wavelength = Wavelength(PredefinedWavelengths.S2A_MSI_08)

    elif BandId == 9:
        s.wavelength = Wavelength(PredefinedWavelengths.S2A_MSI_09)

    elif BandId == 10:
        s.wavelength = Wavelength(PredefinedWavelengths.S2A_MSI_10)

    elif BandId == 11:
        s.wavelength = Wavelength(PredefinedWavelengths.S2A_MSI_11)

    elif BandId == 12:
        s.wavelength = Wavelength(PredefinedWavelengths.S2A_MSI_12)

    elif BandId == 13:
        s.wavelength = Wavelength(PredefinedWavelengths.S2A_MSI_13)

    # 下垫面非均一、朗伯体
    s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromReflectance(-0.1)

    # 运行6s大气模型
    s.run()

    xa = s.outputs.coef_xa
    xb = s.outputs.coef_xb
    xc = s.outputs.coef_xc
    x = s.outputs.values
    return (xa, xb, xc)

def MeanDEM(pointUL, pointDR):
    '''
    计算影像所在区域的平均高程.
    '''
    try:
        DEMIDataSet = gdal.Open("GMTED2km.tif")
    except Exception as e:
        pass

    DEMBand = DEMIDataSet.GetRasterBand(1)
    geotransform = DEMIDataSet.GetGeoTransform()
    # DEM分辨率
    pixelWidth = geotransform[1]
    pixelHight = geotransform[5]

    # DEM起始点：左上角，X：经度，Y：纬度
    originX = geotransform[0]
    originY = geotransform[3]

    # 研究区左上角在DEM矩阵中的位置
    yoffset1 = int((originY - pointUL['lat']) / pixelWidth)
    xoffset1 = int((pointUL['lon'] - originX) / (-pixelHight))

    # 研究区右下角在DEM矩阵中的位置
    yoffset2 = int((originY - pointDR['lat']) / pixelWidth)
    xoffset2 = int((pointDR['lon'] - originX) / (-pixelHight))

    # 研究区矩阵行列数
    xx = xoffset2 - xoffset1
    yy = yoffset2 - yoffset1

    # 读取研究区内的数据，并计算高程
    DEMRasterData = DEMBand.ReadAsArray(xoffset1, yoffset1, xx, yy)

    MeanAltitude = numpy.mean(DEMRasterData)
    return MeanAltitude

def BasicParameters():
    '''
    整理6s大气校正所需的参数
    '''
    global dom
    SixsParameters = dict()
    #太阳天顶角、方位角
    SunAngle = dom.getElementsByTagName('Mean_Sun_Angle')
    SixsParameters["SolarZenithAngle"] = float(SunAngle[0].getElementsByTagName('ZENITH_ANGLE')[0].firstChild.data)
    SixsParameters["SolarAzimuthAngle"] = float(SunAngle[0].getElementsByTagName('AZIMUTH_ANGLE')[0].firstChild.data)

    #卫星天顶角、方位角
    ViewAngles = dom.getElementsByTagName('Mean_Viewing_Incidence_Angle')
    ViewZeniths = dict()
    ViewAzimuths = dict()

    for angle in ViewAngles:
        ViewAngle = int(angle.getAttribute('bandId'))
        ViewZeniths[ViewAngle+1] = float(angle.getElementsByTagName('ZENITH_ANGLE')[0].firstChild.data)
        ViewAzimuths[ViewAngle+1]= float(angle.getElementsByTagName('AZIMUTH_ANGLE')[0].firstChild.data)

    SixsParameters["ViewZenithAngle"] = ViewZeniths
    SixsParameters["ViewAzimuthAngle"] = ViewAzimuths
    # 日期:月、日
    Date = dom.getElementsByTagName('SENSING_TIME')[0].firstChild.data.split('T')[0]
    SixsParameters["ImgMonth"] = int(Date.split('-')[1])
    SixsParameters["ImgDay"] = int(Date.split('-')[2])

    #求影像中心经纬度
    PointULX = int(dom.getElementsByTagName('ULX')[0].firstChild.data)
    PointULY = int(dom.getElementsByTagName('ULY')[0].firstChild.data)

    Imgsizes = dom.getElementsByTagName('Size')

    for Imgsize in Imgsizes:
        Resolution = Imgsize.getAttribute('resolution')
        if Resolution == '10':
            SixsParameters["Nrows"] = int(Imgsize.getElementsByTagName('NROWS')[0].firstChild.data)
            SixsParameters["Ncols"] = int(Imgsize.getElementsByTagName('NCOLS')[0].firstChild.data)
    PointBRX = PointULX + 10*SixsParameters["Ncols"]
    PointBRY = PointULY - 10*SixsParameters["Nrows"]

    # 将投影坐标转为经纬度坐标（具体的投影坐标系由给定数据确定）
    Proj = dom.getElementsByTagName('HORIZONTAL_CS_CODE')[0].firstChild.data
    ProjCode = int(Proj.split(':')[1])

    source = osr.SpatialReference()
    source.ImportFromEPSG(ProjCode)
    target = osr.SpatialReference()
    target.ImportFromEPSG(4326) 
    ct = osr.CoordinateTransformation(source,target)
    CoordsUL,CoordsBR = ct.TransformPoints([(PointULX,PointULY),(PointBRX,PointBRY)])

    ULLon = CoordsUL[0]
    ULLat = CoordsUL[1]
    BRLon = CoordsBR[0]
    BRLat = CoordsBR[1]

    sLongitude = (ULLon+BRLon) / 2
    sLatitude = (ULLat+ULLat) / 2

    #大气模式类型
    if sLatitude > -15 and sLatitude <= 15:
        SixsParameters["AtmosphericProfile"] = 1                                    #Tropical
    elif sLatitude > 15 and sLatitude <= 45:
        if SixsParameters["ImgMonth"] > 4 and SixsParameters["ImgMonth"] <= 9:
            SixsParameters["AtmosphericProfile"] = 2                                #MidlatitudeSummer
        else:
            SixsParameters["AtmosphericProfile"] = 3                                #MidlatitudeWinter
    elif sLatitude > 45 and sLatitude <= 60:
        if SixsParameters["ImgMonth"] > 4 and SixsParameters["ImgMonth"] <= 9:
            SixsParameters["AtmosphericProfile"] = 4                                #SubarctivWinter
        else:
            SixsParameters["AtmosphericProfile"] = 5                                #SubarcticWinter

    pointUL = dict()
    pointDR = dict()
    pointUL["lat"] = ULLat
    pointUL["lon"] = ULLon
    pointDR["lat"] = BRLat
    pointDR["lon"] = BRLon
    SixsParameters["meanDEM"] = (MeanDEM(pointUL, pointDR)) * 0.001

    return SixsParameters

if __name__ == '__main__':

    #输入数据路径
    RootOutName = sys.argv[2]
    RootInputPath = sys.argv[1]

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


            #获得影像头文件 
            MeteData = os.path.join(root,'metadata.xml')
            shutil.copyfile(MeteData,os.path.join(outname,'metedata.xml'))
            dom = xml.dom.minidom.parse(MeteData)
            SixsInputParameter = BasicParameters()

            #选出影像所有波段
            RSbands = glob.glob(os.path.join(root,"B*.tiff"))

            for tifFile in RSbands:
                print(tifFile)
                if os.path.basename(tifFile)=="B8A.tiff":
                    BandId = 9
                elif int(os.path.basename(tifFile)[1:3])<9:
                    BandId = int(os.path.basename(tifFile)[1:3])
                else:
                    BandId = int(os.path.basename(tifFile)[1:3])+1
                print(BandId)
                #捕捉打开数据出错异常
                try:
                    IDataSet = gdal.Open(tifFile)
                except Exception as e:
                    print("文件%S打开失败" % tifFile)
                    LogFile.write('\n'+tifFile+'数据打开失败')
                
                if IDataSet == None:
                    LogFile.write('\n'+tifFile+'数据集读取为空')
                    continue
                else:
                    #获取行列号
                    cols = IDataSet.RasterXSize
                    rows = IDataSet.RasterYSize
                    ImgBand = IDataSet.GetRasterBand(1)
                    ImgRasterData = ImgBand.ReadAsArray(0, 0, cols, rows)
                    if ImgRasterData is None:
                        LogFile.write('\n'+tifFile+'栅格数据为空')
                        continue
                    else:
                        #设置输出文件路径
                        outFilename=os.path.join(outname,os.path.basename(tifFile)[0:3]+'.tiff')

                        #如果文件存在就跳过，进行下一波段操作
                        if os.path.isfile(outFilename):
                            print("%s已经完成" % outFilename)
                            continue
                        else:
                            #表观反射率转换为辐射亮度值
                            RaCalRaster = TOAReflectanceToTOARadiance(BandId)
                            #大气校正
                            a, b, c = AtmosphericCorrection(BandId)
                            y = numpy.where(RaCalRaster!=-9999,a * RaCalRaster - b,-9999)
                            atc = numpy.where(y!=-9999,(y / (1 + y * c))*10000,-9999)
                            
                            driver = gdal.GetDriverByName('GTIFF')
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

    #关闭日志文件
    LogFile.close()


