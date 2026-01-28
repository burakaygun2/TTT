import matplotlib.pyplot as plt 
import numpy as np 
from scipy.optimize import curve_fit

fig, (ax0, ax1) = plt.subplots(1, 2, layout='constrained', figsize=(8, 3)) 

infile = open("data_convection_Titan_80km/statistics.dat", "r") 
lines = infile.readlines() 

time = []
q_top = []
q_bot = []
q_int = []
thickness = []

header = True
for line in lines:
    sline = line.split("\t\t")

    if (header == True):
        header = False
        continue
    else:
        time.append(float(sline[0]))
        q_top.append(float(sline[2])) 
        q_bot.append(float(sline[3])) 
        q_int.append(float(sline[4])) 
        thickness.append(167.0 - float(sline[5]))

ax0.set_xlim(0, time[-1])
ax1.set_xlim(0, time[-1])

ax0.set_ylabel("Heat flux " +r"(mW m$^{-2}$)", labelpad = 5)
ax0.set_xlabel("Time (Myr)", labelpad = 5)
ax0.plot(time, q_top, color="blue", label = "Surface", lw = 2) 
ax0.plot(time, q_bot, color="orange", label = "Ice-water", lw = 2) 
ax0.plot(time, q_int, color="red", label = "Interior", lw = 2) 

ax1.set_ylabel("Ocean thickness (km)", labelpad = 5)
ax1.set_xlabel("Time (Myr)", labelpad = 5)
ax1.plot(time, thickness, color="navy", label = "Ocean thickness", lw = 2) 

ax0.legend(loc="upper right", ncol = 3, fontsize = 8)
plt.savefig("results.png", bbox_inches = "tight", dpi = 250)

# fig, ax0 = plt.subplots(1, 1, layout='constrained', figsize=(4, 4)) 

# for thickness in (80, 120, 160):
#     infile = open("data_convection_Titan_"+str(thickness)+"km/statistics.dat", "r") 
#     lines = infile.readlines() 

#     time = []
#     thickness = []

#     header = True
#     for line in lines:
#         sline = line.split("\t\t")

#         if (header == True):
#             header = False
#             continue
#         else:
#             time.append(float(sline[0]))
#             thickness.append(167.0 - float(sline[5]))

#     # ax0.set_xlim(0, time[-1])
#     ax0.set_ylabel("Ocean thickness (km)", labelpad = 5)
#     ax0.set_xlabel("Time (Myr)", labelpad = 5)
#     ax0.plot(time, thickness, label = r"$D_{ini}=$", lw = 2) 

# ax0.legend(loc="upper right", ncol = 1, fontsize = 8)
# plt.savefig("results.png", bbox_inches = "tight", dpi = 250)
