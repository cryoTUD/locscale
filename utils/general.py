import numpy as np
import math


def round_up_to_even(x):
    return math.ceil(x / 2.) * 2

def round_up_to_odd(x):
    return 1 + round_up_to_even(x)

def true_percent_probability(n):
    x = np.random.uniform(low=0, high=100)
    if x <= n:
        return True
    else:
        return False

