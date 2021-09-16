import numpy as np
import math


def round_up_to_even(x):
    ceil_x = math.ceil(x)
    if ceil_x % 2 == 0:   ## check if it's even, if not return one higher
        return ceil_x
    else:
        return ceil_x+1

def round_up_to_odd(x):
    ceil_x = math.ceil(x)
    if ceil_x % 2 == 0:   ## check if it's even, if so return one higher
        return ceil_x+1
    else:
        return ceil_x

def true_percent_probability(n):
    x = np.random.uniform(low=0, high=100)
    if x <= n:
        return True
    else:
        return False

