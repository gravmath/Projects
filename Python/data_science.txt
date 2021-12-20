### Data Science Module
###
### Copyright - Conrad Schiff (private citizen), June 5, 2021
###
### Version 0.0 - count_switches (no comments)
def count_switches(data):
    import numpy as np
    num_zero          = 0
    num_one           = 0
    index_zero        = []
    index_one         = []
    size_one          = []
    current_index_one = -1
    last              = -1
    #print('i','\t','d','\t','current','\t','last','\t','num_zero','\t','num_one','\t','current_index_one')
    for i,d in zip(range(len(data)),data):
        current = d
        if current != last:
            if d == 0:
                num_zero += 1
                current_index_zero = i
                index_zero.append(current_index_zero)
                if current_index_one != -1:
                    size_one.append(current_index_zero - current_index_one)
            if d == 1:
                num_one  += 1
                current_index_one = i
                index_one.append(current_index_one)
    #    print(i,'\t',d,'\t',current,'\t\t',last,'\t',num_zero,'\t\t',num_one,'\t\t',current_index_one)
        last    = current
    if data[-1] == 1:  #ends with a 1
        size_one.append(i-current_index_one)
    size_one  = np.asarray(size_one)
    index_one = np.asarray(index_one)
    
    return size_one, index_one