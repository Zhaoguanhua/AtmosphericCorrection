#! usr/bin/env python
# -*- coding: utf-8 -*-
# @Author   : Zhaoguanhua
# @DataTime : 2017-10-27 15:08:09

import numpy as np 
# a ="asdfaa"

# print(a[0:3],a[3:5])

b=np.array(([0.5,1.0,0.7,0.45],[-9999,2.3,-0.5,0.9],[2.4,4.5,-9999,0.9],[-9999,3.5,2.7,-9999]))
print(b)
print('\n')
new = np.ones(np.shape(b))
mask = (b!=-9999)
c=b[mask!=0]+1
print(b[mask!=0])
new[mask!=0]=c.astype(float)
print(new)
new[mask==0]=-9999

print(new.astype(float))