import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

#------- 绘图字体设置 ---------#

from matplotlib import rcParams
config = {
    "font.family":'Times New Roman',
    "font.size": 12,
    "mathtext.fontset":'stix',
    "font.serif": ['times'],
}
rcParams.update(config)

#---------- 读入文件 -----------#

data = nc.Dataset('Final/thermocline-20.nc',"r") # [1,70,385,360]
templvl = data['templvl'][:]
depth = data['depth'][:]

#data2 = nc.Dataset('new-new-thermocline-pi0400co21f.nc','r') # [385,360]
#data2 = nc.Dataset('Final/thermocline-10.nc','r') # [385,360]
thermocline = data['thermocline'][:]
plat = data['plat'][:]
plon = data['plon'][:]

#---------- 选择检查的点 -----------#

#[173,0] None
#[174,0] 62.5m

#[134,180] 275m
#[135,180] None

#[260,278] 400m
#[260,279] None

#[356,300]

#----------2022.03.07
#[133,222]
#[130,109]
#[148,152]

#[89,247]

#----------2022.03.22
# xx = 160
# yy = 230

xx = 282
yy = 80

a = thermocline[xx,yy]
templvl = np.array(templvl[0,:,xx,yy])

#-------- 缺省值处理 ---------#

new_templvl = []
new_depth = []
new_a = []

for i in range(0,28) :
    if templvl[i] == -32768:
        break
    new_templvl.append(templvl[i])
    new_depth.append(depth[i])
    new_a.append(a)

#------- 数据展示 --------#

print("")
print("locatation:","%.1f"%plat[xx,yy],"°N, ","%.1f"%plon[xx,yy],"°E")

print("The depth of thermocline is:",a,"m")
print("")

for i in range(18):
    print("--------- ",new_depth[i],"m")
    c = (new_templvl[i] - new_templvl[i+1])/(depth[i+1]-depth[i])
    print("%.4f"%c)

print("------------------")
    
for i in range(18,27):
    print("--------- ",new_depth[i],"m")
    c = (new_templvl[i] - new_templvl[i+1])/(depth[i+1]-depth[i])
    print("%.4f"%c)


#----------- 绘图 ------------#

plt.figure(dpi=120)
ax1 = plt.gca()
ax1.invert_yaxis()

ax1.plot(new_templvl,new_depth)
ax1.plot(new_a)

plt.xlabel('Temperature(deg_C)')
plt.ylabel('Depth(m)')           
ax1.legend(['templvl','thermocline_depth'])

plt.show()

print("")
print('Program Done!')
