def loop_rec(y, number):
   if (number > 1):
      loop_rec( y, number - 1 )
      for i in range(y): 
         print(i, end=' ')        
   else:      
      for i in range(y):
         print(i, end=' ')

loop_rec(4,3)
