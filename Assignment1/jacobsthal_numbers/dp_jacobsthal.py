import sys

memo = {}

def jacobsthal(n):
    if n in memo:
        return memo[n]
    if n == 0:
        return 0
    if n == 1:
        return 1
    memo[n] = jacobsthal(n-1) + 2 * jacobsthal(n-2)
    return memo[n]

sys_n = int(sys.argv[1])
jacobsthal_number = jacobsthal(sys_n)
len_jacobsthal_number = len(str(jacobsthal_number))
print("The ", sys_n, "th Jacobsthal number is: ", jacobsthal_number, "It has ", len_jacobsthal_number, " digits.")