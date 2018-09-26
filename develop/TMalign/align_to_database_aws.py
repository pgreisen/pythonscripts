import sys,os
import shutil
from shutil import copyfile

path='/mys3bucket/databases/'
dir=os.listdir(path)
for i in dir:
        os.chdir(path+'/'+i)
        subdir=os.listdir('./')
        for j in subdir:
                print j
                shutil.copyfile(path+'/'+i+'/'+j, '/home/ec2-user/'+j)
                os.chmod('/home/ec2-user/'+j, 755)
