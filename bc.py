import numpy as np
import numpy.random as random
import os
import sys
from docplex.mp.model import Model



def dropChar(str):
    return str.replace("\n", "")


def readData(path):
    file = open(path ,"r")
    lines = file.readlines() 
    numitem = 0
    cap = 0
    Dalpha = 0
    mu = []
    b = []
    descriptor = ""
    readedlines = 0
    for idx, line in enumerate(lines):
        if readedlines == 0: # descriptor
            descriptor = line
            spltline = line.split(" ")
            numitem = int(dropChar(spltline[-2]))
            cap = float(dropChar(spltline[-3]))
            readedlines += 1
            continue
        if readedlines == 1: # Dalpha
            Dalpha = float(dropChar(line))
            readedlines += 1
            continue
        if readedlines == 2: # mu
            spltline = line.split(" ")
            for c in spltline:
                    if c !=  "\n":
                        mu.append(float(dropChar(c))) 
            readedlines += 1
            continue
        if readedlines == 3: # b
            spltline = line.split(" ")
            for c in spltline:
                if c != "\n":
                    b.append(float(dropChar(c))) 
            readedlines += 1
            break
    #print(numitem, cap, Dalpha, mu, b)
    return numitem, cap, Dalpha, mu, b, descriptor

def best_fit(mus, bs, Dalpha, numitems, cap):
    unpack_items = [i for i in range(numitems)]
    unpack_items.append(-1) # flag
    num_packs = 0
    summu = 0
    sumb = 0
    assignment = [-1 for i in range(numitems)]
    best_fit = -1; 
    skip = -2;
    best_fit_cap_use = cap + 1 #  capacity used
    while unpack_items: 
        item = unpack_items[0]
        unpack_items.pop(0)
        if item == -1: # items in the queue are checked
            if not unpack_items: # every item is assigned
                num_packs += 1
                break
            elif best_fit != -1: # there is an fitting item
                assignment[best_fit] = num_packs
                summu += mus[best_fit]
                sumb += bs[best_fit]
                skip = best_fit
            else: # need new pack
                num_packs += 1
                summu = 0
                sumb = 0
                skip = -2
            unpack_items.append(-1)
            best_fit = -1
            best_fit_cap_use =  cap + 1
            continue;
        if item == skip:
            continue
        lhs = summu + mus[item] + Dalpha * np.sqrt(sumb + bs[item])
        if lhs <= cap: # if fit into the bin
            fit_cap_use = lhs - summu - Dalpha * np.sqrt(sumb) # use of the capacity
            if fit_cap_use < best_fit_cap_use: # update the best fit
                best_fit = item
                best_fit_cap_use = fit_cap_use
        # put the item in checked queue
        unpack_items.append(item)
    sol = [[] for i in range(num_packs)]
    for i in range(numitems):
        sol[assignment[i]].append(i)
    return sol, num_packs
    
def optimize(numitem, cap, Dalpha, lst_mu, lst_b, log_file, result_file,  descriptor, timelimit,  parallelscplex):
    #m = Model()
    m = Model(log_output = log_file)
    sol, numbinpacks = best_fit(lst_mu, lst_b, Dalpha, numitem, cap)
    binpacks = range(numbinpacks)
    # binpack variables
    y = m.binary_var_list(binpacks, name="y_var")
    items = range(numitem)
    x = {(item, binpack): m.binary_var(name="item_var"+str(item)+"."+str(binpack)) for item in items for binpack in binpacks}
    # item times lst_b[item] variables
    xc = {(item, binpack): m.continuous_var(name= str(item)+"."+str(binpack), ub= np.sqrt(lst_b[item])) for item in items for binpack in binpacks}
    # right hand side variables of conic constraints
    z = m.continuous_var_list(binpacks, lb= 0)

    for item in items:
        for binpack in binpacks:
            m.add_constraint(np.sqrt(lst_b[item]) * x[(item, binpack)] <= xc[(item, binpack)])
    for item in items:
        m.add_constraint(m.sum(x[(item, binpack)] for binpack in binpacks) >= 1)
    for binpack in binpacks:
        m.add_constraint(  m.sum(xc[(item, binpack)]**2 for item in items) - z[binpack]**2 <= 0  )
        m.add_constraint(  m.sum(lst_mu[item] * x[(item, binpack)] for item in items) + Dalpha * z[binpack] <= cap*y[binpack])
    m.minimize(m.sum(y))

    # mip start
    warmstart = m.new_solution()
    for binpack in binpacks:
        warmstart.add_var_value(z[binpack], 1)
        for item in sol[binpack]:
            warmstart.add_var_value(x[(item, binpack)], 1)
    m.add_mip_start(warmstart)
    upper_bound = len(binpacks)

    # parameter
    m.parameters.clocktype = 1
    m.set_time_limit(timelimit)
    m.parameters.threads.set(0 if parallelscplex  else 1)
    m.parameters.mip.tolerances.absmipgap.set(1.0)
    #m.parameters.mip.strategy.search.set(2)
    s = m.solve()
    if s:
        f = result_file
        m.report()
        status = m.solve_status 
        soltime =  m.solve_details.time
        best_bound =  m.solve_details.best_bound
        obj = s.get_objective_value()
        nb_nodes =  m.solve_details.nb_nodes_processed
        relgap = round(abs(best_bound - obj) / (best_bound + 1e-6) , 4) *100
        relgap_round = round(abs(np.ceil(best_bound) - obj) / (np.ceil(best_bound) + 1e-6) , 4) * 100
        f.write(descriptor +  '\n')
        f.write('status: ' + str(status) + '\n')
        f.write('soltime: ' + str(soltime) + '\n')
        f.write('nodes: ' + str(nb_nodes) + '\n')
        f.write('relgap: ' + str(relgap) + '\n')
        f.write('obj: ' + str(obj) + '\n')
        f.write('best_bound: ' + str(best_bound) + '\n')
        f.write('relgap_round: ' + str(relgap_round) + '\n')
        f.write('first_solution: ' + str(upper_bound) + '\n')
        f.close()
    else:
        print("* model has no solution")
    log_file.close()


def main(argv):
    path = os.getcwd()
    timelimit = float(argv[1])
    parallelscplex = int(argv[2]) == 1
    dataname= argv[3]
    data_rpath= argv[4]
    log_rpath= argv[5]
    result_rpath= argv[6]
    data_path = os.getcwd() +  "/"+ data_rpath +   "/"+dataname 
    log_path = os.getcwd() +  "/"+ log_rpath +  "/"+dataname + "_bcsocp.log"
    result_path = os.getcwd() +  "/"+ result_rpath +   "/"+dataname + "_bcsocp.result"
    numitem, cap, Dalpha, mu, b,   descriptor = readData(data_path)
    log_file = open(log_path, "w")
    result_file =  open(result_path, "w")
    optimize(numitem, cap, Dalpha, mu, b, log_file, result_file,  descriptor, timelimit,  parallelscplex)
 
if __name__ == "__main__":
    main(sys.argv)