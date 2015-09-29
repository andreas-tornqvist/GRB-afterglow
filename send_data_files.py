import os
for i in range(3):
    for j in range(3):
        for k in range(3):
            os.system('scp -rP 1047 ~/Data_files/fullSED_%d%d%dfit/ $toSol:Data_files/'%(i,j,k))
