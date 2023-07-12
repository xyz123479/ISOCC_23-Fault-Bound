from multiprocessing import Pool
import sys
import os
import time

oecc = [0,1,2,4,5,6,7,8,9,10,11] # 0 : 1 x 16 (16b bound), 1 : 2 x 8 (16b bound), 2 : 4 x 4 (16b bound), 3 : 1 x 32 (32b bound), 4 : 2 x 16 (32b bound), 5 : 4 x 8 (32b bound), 6 : Un_bound_1_16, 7 : Un_bound_2_8, 8 : Un_bound_4_4, 9 : Un_bound_1_32, 10 : Un_bound_2_16, 11 : Un_bound_4_8
fault = [5,6,7,8,12,13] # SE_SE=0, SE_MBBE=1, SE_SW=2, SE_SP=3, SE_CHIPKILL=4, MBBE_MBBE=5, MBBE_SW=6, MBBE_SP=7, MBBE_CHIPKILL=8, SW_SW=9, SW_SP=10, SW_CHIPKILL=11, SP_SP=12, SP_CHIPKILL=13, CHIPKILL_CHIPKILL=14
recc = [0,1,2,3,4,5] # 0 : 4 x 2 SSC, 1 : 4 X 4 SSC, 2 : 2 x 4 DSC, 3 : 2 x 8 DSC, 4 : 1 x 8 QSC, 5 : 1 x 16 QSC, 6 : No-Rank level ECC

# SE: per-chip Single bit Error
# MBBE: per-chip Muli bit Bounded Error
# SW: per-chip Single Word Error (all 4 bits)
# SP: per-chip Single Pin Error (More than 2 bits)
# CHIPKILL : Chip error (All random)

for oecc_param in oecc:
    for fault_param in fault:
        for recc_param in recc:
            os.system("./Fault_sim_start {0:d} {1:d} {2:d} &".format(oecc_param, fault_param, recc_param))
