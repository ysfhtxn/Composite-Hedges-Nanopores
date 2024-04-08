import numpy as np
from itertools import product

a = np.zeros([5,5], dtype='str')
a[:] = 'a'
print(a)

position_list = list(product(range(5), range(5)))
print(position_list )