#!/usr/bin/env python3
# fibonacci.py

"""Question 6"""

import sys

# pylint: disable=redefined-outer-name
# pylint: disable=consider-using-enumerate
# pylint: disable=invalid-name
# pylint: disable=unused-variable


def population(n, k):
    """ This function takes a day and a reproduction
    rate and returns the population
    size at day n
    :param n: day
    :param k: reproduction rate
    :return: current
    """
    # set values for yesterday and day before that
    prev_1, prev_2 = 1, 0

    if n == 0:
        return prev_2
    if n == 1:
        return prev_1

    for i in range(n - 1):
        current = prev_1 + k * prev_2
        # reassign values to yesterday, day before, and today
        prev_1, prev_2 = current, prev_1

    return current


if __name__ == "__main__":
    # Check to make sure there are at least two arguments

    arg_count = len(sys.argv) - 1
    if arg_count < 2:
        raise Exception("This script requires 2 arguments: day and "
                        "reproduction rate")

    n = int(sys.argv[1])
    k = int(sys.argv[2])

    print(population(n, k))
