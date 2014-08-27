#!/usr/bin/env python
import os,shutil

'''

If file is present it moves the directory to a destination


'''


def main():

    path = './'
    
    dst = 'Destination'

    tricker_file = '0score.sc'

    os.mkdir(dst)

    dirs = os.listdir(path)
    
    for dir in dirs:
    
        if os.path.isdir(dir):

            os.chdir(dir)
            
            copy_dir = False
            
            files = os.listdir(path)
            
            for fl in files:
                
                if(fl == tricker_file):
                    
                    copy_dir = True
                    print dir
            os.chdir('../')

            #if(copy_dir):   
            #    shutil.copytree(dir,dst+'/'+dir)

      
if __name__ == "__main__":
    main()
