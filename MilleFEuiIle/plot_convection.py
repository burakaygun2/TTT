import matplotlib.pyplot as plt 
import numpy as np 
from scipy.optimize import curve_fit

from matplotlib.patches import Rectangle


fig, ax0 = plt.subplots(1, 1, layout='constrained', figsize=(6, 4)) 

thickness_list = [20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 167]
viscosity_list = [13, 14, 15, 16]

ax0.set_xlim(10, 175)
ax0.set_ylim(12.5, 16.99)

ax0.set_yticks([13, 14, 15, 16], [r"$10^{13}$", r"$10^{14}$", r"$10^{15}$", r"$10^{16}$"])
ax0.set_xticks([20, 40, 60, 80, 100, 120, 140, 167])

data_x_conv = []
data_y_conv = []
data_c_conv = []

data_x_cond = []
data_y_cond = []
data_c_cond = []

diss_max = 0
diss_min = 10   
for thickness in thickness_list:
    for viscosity in viscosity_list:
        try:
            if (viscosity == 16):
                filename = "/nfs90/kihoulou/TTT/data_Titan_" + str(thickness) + "km_" + str(viscosity) + "Pas/statistics.dat"
            elif (viscosity == 13):
                filename = "/nfsy2/kihoulou/TTT/data_Titan_" + str(thickness) + "km_" + str(viscosity) + "Pas/statistics.dat"
            else:
                filename = "/nfs00/kihoulou/TTT/data_Titan_" + str(thickness) + "km_" + str(viscosity) + "Pas/statistics.dat"
            
            if (thickness == 30 or thickness == 20):
                filename = "/nfs00/kihoulou/TTT/data_Titan_" + str(thickness) + "km_" + str(viscosity) + "Pas/statistics.dat"

            infile = open(filename, "r") 
            lines = infile.readlines() 

            time = []
            q_top = []
            vrms = []

            header = True
            for line in lines:
                sline = line.split("\t\t")

                if (header == True):
                    header = False
                    continue
                else:
                    time.append(float(sline[0]))
                    q_top.append(float(sline[2])/1e3)
                    vrms.append(float(sline[4]))

            qtop = []
            for i in range(0, len(time)-1):
                if (time[i] > time[-1] - 1):
                    qtop.append(q_top[i])

            qtop_aver = sum(qtop) / len(qtop)

            if (vrms[0] < vrms[-1]):
                data_x_conv.append(thickness)
                data_y_conv.append(viscosity)
                data_c_conv.append(qtop_aver*4.0*np.pi*2570e3**2/1e12)

            else:
                data_x_cond.append(thickness)
                data_y_cond.append(viscosity)
                data_c_cond.append(qtop_aver*4.0*np.pi*2570e3**2/1e12)

            if (qtop_aver*4.0*np.pi*2570e3**2/1e12 > diss_max):
                diss_max =  qtop_aver*4.0*np.pi*2570e3**2/1e12

            if (qtop_aver*4.0*np.pi*2570e3**2/1e12 < diss_min):
                diss_min =  qtop_aver*4.0*np.pi*2570e3**2/1e12

            ax0.text(thickness - 3,viscosity + 0.15, round(qtop_aver*4.0*np.pi*2570e3**2/1e12, 1), fontsize=8)
        except:
            continue

ax0.set_ylabel("Basal ice shell viscosity (Pa s)", labelpad = 5)
ax0.set_xlabel("Ice shell thickness (km)", labelpad = 5)

def ice_thickness(x):
    return 167 - x

def ocean_thickness(x):
    return 167 - x

secax = ax0.secondary_xaxis('top', functions=(ice_thickness, ocean_thickness))
secax.set_xlabel('Ocean thickness (km)', labelpad = 5)
secax.set_xticks([0, 20, 40, 60, 80, 100, 120, 140, 167])


ax0.scatter(data_x_cond, data_y_cond, c = data_c_cond, cmap = "coolwarm", s = 40, marker = "o", vmin = round(diss_min, 1), vmax = round(diss_max, 1))
cc = ax0.scatter(data_x_conv, data_y_conv, c = data_c_conv, cmap = "coolwarm", s = 40, marker = "v", vmin = round(diss_min, 1), vmax =round(diss_max, 1)) 


ax0.scatter(0, 0, c ="black", s = 40, marker = "o", label = "Conduction") 
ax0.scatter(0, 0, c ="black", s = 40, marker = "v", label = "Convection") 
thickness_list = np.linspace(10, 180)

constant=81e3**3/1e15
data_x_fit = []
for i in range(0, len(thickness_list)):
    data_x_fit.append(np.log10((thickness_list[i]*1e3)**3/constant))

ax0.plot(thickness_list, data_x_fit, color="gray", lw = 1, linestyle="dashed", dashes = [5,5], zorder = 0) 

cbar= plt.colorbar(cc, label="Surface heat flux (TW)", ticks=np.linspace(round(diss_min, 1), round(diss_max, 1), 4))

ax0.legend(loc="upper right", ncol = 2, fontsize = 8)
ax0.text(15, 14.5, "Conduction", c="grey")
ax0.text(70, 14.5, "Convection", c="grey")

rect = Rectangle([35, 12.75], 140, 0.75, color="white", zorder = 20, alpha = 0.9)
ax0.add_patch(rect)

plt.savefig("Fig1.png", bbox_inches = "tight", dpi = 250)
plt.savefig("Fig1.pdf", bbox_inches = "tight")
