import numpy as np
import itertools
from scipy.linalg import expm

x = 'ACGTA'
y = 'AAATA'
T = 0.02
phi = 10.0

def base_to_numeric(sequence):
    base_mapping = {'A':0,'G':1,'C':2,'T':3}
    return np.array([base_mapping[i] for i in sequence])

def numeric_to_base(sequence):
    base_mapping = {0:'A',1:'G',2:'C',3:'T'}
    return ''.join([base_mapping[i] for i in sequence])

x = base_to_numeric(x)
y = base_to_numeric(y)


def rate(x_new,x):
    if (x[0] == 2) & (x[1] == 1):
        return 1/3 * phi
    elif (x[1] == 2) & (x[2] == 1):
        return 1/3 * phi
    else:
        return 1/3

bases = [0, 1, 2, 3]
combinations = list(itertools.product(bases, repeat=3))
seqs = [np.array(combination) for combination in combinations]

def row(i):
    row_rates = []
    for j in seqs:
        if (np.sum(i != j) != 1):
            row_rates.append(0)
        else:
            diff = [k for k in range(len(i)) if i[k] != j[k]]
            x_new = j[diff[0]]
            if diff[0] == 0:
                x = [0,i[0],i[1]]
            elif diff[0] == 2:
                x = [i[1],i[2],0]
            else:
                x = i
            row_rates.append(rate(x_new,x))
    return row_rates

x_trim = x[1:-1]
y_trim = y[1:-1]

x_index = [k for k in range(len(seqs)) if np.array_equal(seqs[k],x_trim)][0]
y_index = [k for k in range(len(seqs)) if np.array_equal(seqs[k],y_trim)][0]

Q = np.zeros((64, 64))
for i in range(len(seqs)):
    rates = row(seqs[i])
    Q[i] = rates
    Q[i,i] = -np.sum(rates)

def forward_simulation(Q,T):
    times = [0]
    states = [x_index]
    positive = True
    while positive:
        new_time = np.random.exponential(-1/Q[states[-1],states[-1]])
        times.append(times[-1] + new_time)
        if (T - times[-1]) > 0:
            possible_states = np.arange(0,64)[np.arange(0,64) != states[-1]]
            prob = Q[states[-1],possible_states] / -Q[states[-1],states[-1]]
            states.append(np.random.choice(possible_states,size=1,replace=False,p=prob)[0])
        else:
            states.append(states[-1])
            positive = False
    return times,states

def accept_reject(Q,T):
    sample_num = 8192
    samples = []
    while sample_num > 0:
        fs = forward_simulation(Q,T)
        if fs[1][-1] == y_index:
            sample_num = sample_num - 1
            samples.append(fs)
    return samples

all_samples = accept_reject(Q,T)

with open('example_output/true_number_of_mutations.txt', 'w') as f:
    for i in range(len(all_samples)):
        times, states = all_samples[i]
        f.write(str(len(times) - 2))
        f.write('\n')


with open('example_output/true_first_jump_time.txt', 'w') as f:
    for i in range(len(all_samples)):
        times, states = all_samples[i]
        f.write(str(times[1]))
        f.write('\n')