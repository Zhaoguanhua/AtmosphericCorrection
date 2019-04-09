#! usr/bin/env python
# -*- coding:utf-8 -*-
# created by zhaoguanhua 2017/10/23

from osgeo import gdal
from osgeo import osr
import numpy as np

def getSRSPair(dataset):
    '''
    获得给定数据的投影参考系和地理参考系
    :param dataset: GDAL地理数据
    :return: 投影参考系和地理参考系
    '''
    prosrs = osr.SpatialReference()

    prosrs.ImportFromWkt(dataset.GetProjection())
    geosrs = prosrs.CloneGeogCS()
    return prosrs, geosrs

def geo2lonlat(dataset, x, y):
    '''
    将投影坐标转为经纬度坐标（具体的投影坐标系由给定数据确定）
    :param dataset: GDAL地理数据
    :param x: 投影坐标x
    :param y: 投影坐标y
    :return: 投影坐标(x, y)对应的经纬度坐标(lon, lat)
    '''
    prosrs, geosrs = getSRSPair(dataset)
    ct = osr.CoordinateTransformation(prosrs, geosrs)
    coords = ct.TransformPoint(x, y)
    return coords[:2]

def lonlat2geo(dataset, lon, lat):
    '''
    将经纬度坐标转为投影坐标（具体的投影坐标系由给定数据确定）
    :param dataset: GDAL地理数据
    :param lon: 地理坐标lon经度
    :param lat: 地理坐标lat纬度
    :return: 经纬度坐标(lon, lat)对应的投影坐标
    '''
    prosrs, geosrs = getSRSPair(dataset)
    ct = osr.CoordinateTransformation(geosrs, prosrs)
    coords = ct.TransformPoint(lon, lat)
    return coords[:2]

def imagexy2geo(dataset, row, col):
    '''
    根据GDAL的六参数模型将影像图上坐标（行列号）转为投影坐标或地理坐标（根据具体数据的坐标系统转换）
    :param dataset: GDAL地理数据
    :param row: 像素的行号
    :param col: 像素的列号
    :return: 行列号(row, col)对应的投影坐标或地理坐标(x, y)
    '''
    trans = dataset.GetGeoTransform()
    px = trans[0] + row * trans[1] + col * trans[2]
    py = trans[3] + row * trans[4] + col * trans[5]
    return px, py

def geo2imagexy(dataset, x, y):
    '''
    根据GDAL的六 参数模型将给定的投影或地理坐标转为影像图上坐标（行列号）
    :param dataset: GDAL地理数据
    :param x: 投影或地理坐标x
    :param y: 投影或地理坐标y
    :return: 影坐标或地理坐标(x, y)对应的影像图上行列号(row, col)
    '''
    trans = dataset.GetGeoTransform()
    a = np.array([[trans[1], trans[2]], [trans[4], trans[5]]])
    b = np.array([x - trans[0], y - trans[3]])
    return np.linalg.solve(a, b)  # 使用numpy的linalg.solve进行二元一次方程的求解

if __name__ == '__main__':
    gdal.AllRegister()
    dataset = gdal.Open('/Users/zhaoguanhua/2017/Project/Data/AtmosphericCorrection/sentinel/15SUC/tiles/15/S/UC/2017/1/7/0/B02.tiff')
    print('数据投影：')
    print(dataset.GetProjection())
    print('数据的大小（行，列）：')
    print('(%s %s)' % (dataset.RasterYSize, dataset.RasterXSize))

    x = 300000
    y = 4300020
    lon = 122.47242
    lat = 52.51778
    row = 2399
    col = 3751

    print('投影坐标 -> 经纬度：')
    coords = geo2lonlat(dataset, x, y)
    print('(%s, %s)->(%s, %s)' % (x, y, coords[0], coords[1]))
    print('经纬度 -> 投影坐标：')
    coords = lonlat2geo(dataset, lon, lat)
    print('(%s, %s)->(%s, %s)' % (lon, lat, coords[0], coords[1]))

    print('图上坐标 -> 投影坐标：')
    coords = imagexy2geo(dataset, row, col)
    print('(%s, %s)->(%s, %s)' % (row, col, coords[0], coords[1]))
    print('投影坐标 -> 图上坐标：')
    coords = geo2imagexy(dataset, x, y)
    print('(%s, %s)->(%s, %s)' % (x, y, coords[0], coords[1]))




