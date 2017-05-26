import os
import threading
import time
import glob

def find_orb():
    os.system("./fo mpc.txt")

path= 'MPC'
dirs = os.listdir("MPC")
mpcArray= []
imageArray= []
i=0

for filename in glob.glob(os.path.join(path, '*.txt')):
    r= open(filename, "r")
    for line in r.readlines():
        mpcArray.append(line.strip())
        i=i+1
    r.close()

f= open("mpc.txt", "w")
g= open("image.txt", "w")
k= 0
size=0

for j in range(i):
    f.write("     "+mpcArray[j]+'\n')
    size=size+1

f.close()
g.close()
print(size)
t = threading.Thread(target= find_orb)
t.daemon= True
t.start()

time.sleep(size/4)
os.system("killall -9 fo")
os.remove('MPCORB.DAT')
os.rename('mpc_fmt.txt','MPCORB.DAT')
os.system('cp -u MPCORB.DAT /mnt/c/Astrometrica/Catalogs')
