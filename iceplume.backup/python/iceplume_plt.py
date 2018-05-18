import numpy as np
import csv
# read from output file
fh = open('../iceplume_test_output.txt', 'r')

zw = np.linspace(-40., 0., 41)
data = csv.reader(fh, delimiter=' ', skipinitialspace=True)

rad = []
vel = []
temp = []
salt = []
area = []
melt = []
for row in data:
    rad.append(float(row[0]))
    vel.append(float(row[1]))
    temp.append(float(row[2]))
    salt.append(float(row[3]))
    area.append(float(row[4]))
    melt.append(float(row[5]))

# fh.close()
# for i in range(41):
#     fh.read('%f, %f, %f\n' % (zw[i], saltr[i], temp[i]))
# fh.close()
