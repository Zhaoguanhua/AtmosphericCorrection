#! usr/bin/env python
# -*- coding: utf-8 -*-
# @Author   : Zhaoguanhua
# @DateTime : 2018-01-05 11:05:11

from Py6S import *
import numpy as np

def testpy6s():
    # 6S模型
    s = SixS()

    # 传感器类型 自定义
    s.geometry = Geometry.User()
    s.geometry.solar_z = 36
    s.geometry.solar_a = 153
    s.geometry.view_z = 0
    s.geometry.view_a = 0

    # 日期
    s.geometry.month = 3
    s.geometry.day = 27

    #print(s.geometry)

    s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeWinter)

    #print(s.atmos_profile)
    # 气溶胶类型大陆
    s.aero_profile = AtmosProfile.PredefinedType(AeroProfile.Continental)

    # 下垫面类型
    s.ground_reflectance = GroundReflectance.HomogeneousLambertian(0.36)

    # 550nm气溶胶光学厚度,对应能见度为40km
    s.aot550 = 0.14497
    # s.aot550 = MeanAot(Dateparm[1])

    # 研究区海拔、卫星传感器轨道高度
    s.altitudes = Altitudes()
    s.altitudes.set_target_custom_altitude(0.059)
    s.altitudes.set_sensor_satellite_level()

    # 校正波段（根据波段名称）
    s.wavelength = Wavelength(0.450,0.520,[0.26665045403732035, 0.3898905990356667, 0.4365734433632099, 0.46350262095778305, 0.49204568689798206, 0.5176199163843476, 0.540561115963024, 0.5682131273380822, 0.6002756993408136, 0.6221870897658232, 0.6305517292683029, 0.6391171618716124, 0.6629074375761411, 0.6941437481584642, 0.7198332777418507, 0.7410571188969691, 0.7656783699869292, 0.7967537861568459, 0.8253006271461241, 0.8405207257502543, 0.8417219636907287, 0.8447305755272894, 0.8699114324074595, 0.9136765723097032, 0.9508133004267264, 0.9819957569853817, 0.9946839532212262, 0.8918374093071325, 0.5940393710297492])
                        
    print(s.wavelength)


    # s.atmos_corr = AtmosCorr.AtmosCorrBRDFFromRadiance()
    # s.atmos_corr = AtmosCorr.AtmosCorrBRDFFromReflectance(15)
    # s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromRadiance(0.1)
    s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromReflectance(-0.5)

    print(s.atmos_corr)

    # 运行6s大气模型
    s.run()
    xa = s.outputs.coef_xa
    xb = s.outputs.coef_xb
    xc = s.outputs.coef_xc
    x = s.outputs.values
    print(xa,xb,xc)
    print(x)
    return (xa, xb, xc)

def testL8py6s():
    # 6S模型
    s = SixS()

    # 传感器类型 自定义
    s.geometry = Geometry.User()
    s.geometry.solar_z = 46
    s.geometry.solar_a = 155
    s.geometry.view_z = 0
    s.geometry.view_a = 0

    # 日期
    s.geometry.month = 11
    s.geometry.day = 1

    #print(s.geometry)

    s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeWinter)

    #print(s.atmos_profile)
    # 气溶胶类型大陆
    s.aero_profile = AtmosProfile.PredefinedType(AeroProfile.Continental)

    # 下垫面类型
    s.ground_reflectance = GroundReflectance.HomogeneousLambertian(0.36)

    # 550nm气溶胶光学厚度,对应能见度为40km
    s.aot550 = 0.14497

    # 研究区海拔、卫星传感器轨道高度
    s.altitudes = Altitudes()
    s.altitudes.set_target_custom_altitude(0.112)
    s.altitudes.set_sensor_satellite_level()

    # 校正波段（根据波段名称）
    s.wavelength = Wavelength(0.436, 0.526,
                      [1.00000000e-05, 1.79000000e-04, 4.55000000e-04,
                                1.63350000e-03, 6.86900000e-03, 4.28880000e-02,
                                2.71370000e-01, 7.90740500e-01, 9.03034000e-01,
                                9.04677500e-01, 8.89667000e-01, 8.79232000e-01,
                                8.79688000e-01, 8.89796500e-01, 8.48533000e-01,
                                8.36270500e-01, 8.68497000e-01, 9.11461500e-01,
                                9.31726000e-01, 9.54896500e-01, 9.56424000e-01,
                                9.83834000e-01, 9.89469000e-01, 9.68066500e-01,
                                9.88729000e-01, 9.61057500e-01, 9.66125000e-01,
                                9.82077000e-01, 9.63135000e-01, 9.98249000e-01,
                                8.44893000e-01, 1.19533500e-01, 5.32800000e-03,
                                1.32850000e-03, 5.16000000e-04, 1.17000000e-04,
                                2.30000000e-05])
    # s.wavelength = Wavelength(PredefinedWavelengths.LANDSAT_OLI_B2)
    print(s.wavelength)


    # s.atmos_corr = AtmosCorr.AtmosCorrBRDFFromRadiance()
    # s.atmos_corr = AtmosCorr.AtmosCorrBRDFFromReflectance(15)
    # s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromRadiance(0.1)
    s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromReflectance(-0.1)

    print(s.atmos_corr)

    # 运行6s大气模型
    s.run()
    xa = s.outputs.coef_xa
    xb = s.outputs.coef_xb
    xc = s.outputs.coef_xc
    x = s.outputs.values
    print(xa,xb,xc)
    print(x)
    return (xa, xb, xc)

if __name__ == '__main__':
    # testpy6s()
    testL8py6s()

