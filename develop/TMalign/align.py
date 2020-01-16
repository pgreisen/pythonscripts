from pymol import cmd
#run /Users/pgreisen/pymolscripts/tmalign.py
all_objs = cmd.get_names("objects",enabled_only=0)
for i in cmd.get_object_list():
    print "Aligning the following: ",i
    tmalign( i, '3vte' ) 
    cmd.center( i ,animate=-1)
