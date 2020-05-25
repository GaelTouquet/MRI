def make_fibonacci(length):
     if length == 1:
          return [1]
     fibo = [1,2]
     for i in range(length-2):
          fibo.append(fibo[-2] + fibo[-1])
     return fibo

def find_best_number_of_points(Ntarget, fibo):
    i = 0
    previous_number = 0
    new_number = 0
    while (i * fibo) < Ntarget:
        i += 1
        previous_number = new_number
        new_number = i * fibo
    if abs(previous_number - Ntarget) < abs(new_number - Ntarget):
        return i - 1 
    else:
        return i

def find_best_fibo(Ntarget, n_readouts_per_interleaf, fibonacci):
    previous_number = 0
    new_number = 0
    for i in range(len(fibonacci)):
        previous_number = new_number
        new_number = fibonacci[i] * n_readouts_per_interleaf
        if new_number > Ntarget:
            break
    if abs(previous_number - Ntarget) < abs(new_number - Ntarget):
        return fibonacci[i-1]
    else:
        return fibonacci[i]