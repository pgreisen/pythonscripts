'''
Get missing values from hash-table

'''

A = [1,2,3]
B = [3,4]
C = [1,3,4]
D = [4,5]

centres = []
missing_value = {}

centres.append(A)
centres.append(B)
centres.append(C)
centres.append(D)

dummy = 0

for i in centres:
    
    for j in i:
        missing_value[j] = dummy

    dummy += 1

print missing_value

print "A: ",A
print "B: ",B
print "C: ",C
print "D: ",D

# Fill out missing values
for center in centres:
    for key in missing_value.keys():
        if key not in center:
            center.append(key)

print "#################################"
print "A: ",A
print "B: ",B
print "C: ",C
print "D: ",D
