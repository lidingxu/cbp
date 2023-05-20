import numpy as np
import numpy.random as random
from scipy.stats import norm
from scipy.stats import truncnorm
import os
import copy
import time

path = os.getcwd()

from pathlib import Path
import os

benchmarks = ["TrainMedium", "TrainBig"]
numitems =  { "TrainMedium": 500, "TrainBig": 900}

caps = [72]                        
alphas = [0.6,0.7,0.8,0.9,0.95,0.99]
cases =  ["g", "h", "d"]
seeds = []
num_seed = 2




def get_truncated_normal(mean=0, sd=1, low=0, upp=10):
    trunc =  truncnorm(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)
    mean = trunc.moment(1) 
    var = trunc.moment(2)
    return mean, np.sqrt(var)

def sample_item():
    dict_param = {}
    # A_lb uniform in [0.3, 0.6]
    dict_param["A_lb"] = random.uniform(0.3, 0.6)
    # A_ub uniform in [0.7, 1.0]
    dict_param["A_ub"] = random.uniform(0.7, 1.0)
    # mu uniform in [0.1, 0.5]
    mu_ = random.uniform(0.1, 0.5)
    # sigma uniform in [0.1, 0.5]   
    sigma_ = random.uniform(0.1, 0.5)
    mu, sigma = get_truncated_normal(mu_, sigma_, dict_param["A_lb"], dict_param["A_ub"])
    dict_param["mu"] = mu
    dict_param["sigma"] = sigma

    return dict_param


def gaussian_case(lst_mu, lst_param,  alpha):
    lst_b = []
    for mu, dict_param in zip(lst_mu, lst_param):
        lst_b.append((dict_param["sigma"] * mu)**2)
    Dalpha = norm.ppf(alpha)


    return lst_b, Dalpha


def heoffding_inequality(lst_mu, lst_param,  alpha):
    lst_b = []
    for mu, dict_param in zip(lst_mu, lst_param):
        lst_b.append((dict_param["A_ub"]  * mu - dict_param["A_lb"]  * mu )**2)
    Dalpha = np.sqrt(-np.log(1-alpha)/2)


    return lst_b, Dalpha

def distribute_ro(lst_mu, lst_param,  alpha):
    lst_b = []
    for mu, dict_param in zip(lst_mu, lst_param):
        lst_b.append((dict_param["sigma"]  * mu)**2)
    Dalpha = np.sqrt(alpha/ (1 -alpha))


    return lst_b, Dalpha

def find_ind(prob_distr, ele):
    ind = 0
    for i, prob in enumerate(prob_distr):
        if ele >= prob[0] and ele <prob[1]:
            ind = min(max(i, 0), len(prob_distr) - 1)
            break
    return ind

def sample_itm(numitem):
    prob_distr_ = [36.3/100, 13.8 / 100, 21.3/100, 23.1/100, 3.5/100, 1.9/100]
    prob_distr = [36.3/100, 13.8 / 100, 21.3/100, 23.1/100, 3.5/100, 1.9/100]
    dens_ = 0
    for i, dens in enumerate(prob_distr):
        dens_ += prob_distr[i] 
        prob_distr[i] = (dens_ - prob_distr[i], dens_  + (0.1 if i == len(prob_distr) -1 else 0))
    sizesitm = [1, 2 ,4, 8, 16, 32]
    numsize = len(sizesitm)
    print(prob_distr)
    lmu = []
    for i in range(numitem):
        rd_sample = random.uniform()
        ind = find_ind(prob_distr, rd_sample)
        inter = prob_distr[ind]
        dp = rd_sample - inter[0]
        if ind != len(sizesitm) -1:
            dsize = sizesitm[min(ind + 1, numsize -1)] - sizesitm[ind] 
        else:
            dsize = 1
        size = dp / prob_distr_[ind] * dsize + sizesitm[ind] 
        lmu.append(size)
    return lmu

def gen_instance(numitem, cap, alpha, case):
    if case == "g":
        gen = gaussian_case
    elif case == "h":
        gen = heoffding_inequality
    elif case == "d":
        gen = distribute_ro
    lst_param = []
    for _ in range(numitem):
        lst_param.append(sample_item())
    lst_mu = sample_itm(numitem)
    lst_b, Dalpha = gen(lst_mu, lst_param, alpha)
    for i in range(numitem):
        #print( lst_param[i]["mu"] )
        lst_mu[i] = lst_mu[i]* lst_param[i]["mu"] 
        ratio = max((lst_mu[i] + Dalpha * np.sqrt(lst_b[i]) + 5) / cap,  1)
        lst_mu[i] /= ratio
        lst_b[i] /= ratio * ratio
    return lst_mu, lst_b, Dalpha


def first_fit(lst_mu, lst_b, Dalpha, numitem, cap):
    unfixed = [i for i in range(numitem)]
    numbins = 0
    summu = 0
    sumb = 0
    while len(unfixed) > 0:
        unfixed_ = copy.deepcopy(unfixed)
        for item in unfixed_:
            lhs = summu + lst_mu[item] + Dalpha * np.sqrt(sumb + lst_b[item])
            if lhs <= cap:
                summu += lst_mu[item]
                sumb += lst_b[item]
                unfixed.remove(item)
        numbins += 1
        summu = 0
        sumb = 0

    return numbins

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


def test_bin_pack_cplex(lst_mu, lst_b, Dalpha, numitem, cap, file_name, descriptor):
    flog = open(path+ "/logs/"+file_name + "_bin.log", "w")
    m = Model(log_output = flog)
    sol, numbinpacks = best_fit(lst_mu, lst_b, Dalpha, numitem, cap)
    binpacks = range(numbinpacks)
    # binpack variables
    y = m.binary_var_list(binpacks)
    items = range(numitem)
    # item variables
    x = {(item, binpack): m.binary_var() for item in items for binpack in binpacks}
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

    # parameter
    m.parameters.clocktype = 1
    m.set_time_limit(3600)
    s = m.solve()
    if s:
        f = open(path + "/reports/" + file_name + ".report","x")
        m.report()
        status = m.solve_status 
        soltime =  m.solve_details.time
        best_bound =  m.solve_details.best_bound
        obj = s.get_objective_value()
        nb_nodes =  m.solve_details.nb_nodes_processed
        relgap = round(abs(best_bound - obj) / (best_bound) , 4) *100
        relgap_round = round(abs(np.ceil(best_bound) - obj) / (np.ceil(best_bound) ) , 4) * 100
        f.write(descriptor +  '\n')
        f.write(descriptor +  '\n')
        f.write('status: ' + str(status) + '\n')
        f.write('soltime: ' + str(soltime) + '\n')
        f.write('nodes: ' + str(nb_nodes) + '\n')
        f.write('relgap: ' + str(relgap) + '\n')
        f.write('obj: ' + str(obj) + '\n')
        f.write('best_bound: ' + str(best_bound) + '\n')
        f.write('relgap_round: ' + str(relgap_round) + '\n')
        f.close()
    else:
        print("* model has no solution")

def checkpackable(mus, bs, Dalpha, numitems, cap):
    for i in range(numitems):
        if mus[i] + Dalpha * np.sqrt(bs[i]) > cap:
            print(mus[i], Dalpha , np.sqrt(bs[i]))
            return False
    return True

for _ in range(num_seed):
    seeds.append(random.randint(1000))

data_path = os.getcwd() + "/train_data"

cap = 72
for benchmark in benchmarks:
    int_bench = 0
    for char in benchmark:
        int_bench = ord(char)
    numitem = numitems[benchmark]
    for seed in seeds:
        np.random.seed(seed + int_bench)
        for case in cases:
            for alpha in alphas:
                packable = False
                while not packable: 
                    lst_mu, lst_b, Dalpha = gen_instance(numitem, cap, alpha, case)
                    packable = checkpackable(lst_mu, lst_b, Dalpha, numitem, cap)
                descriptor =  str(case) + " " +  str(cap)  +  " "  + str(numitem) + " " + str(alpha)
                file_name = str(case) + "_" +  str(cap)  +  "_"  + str(numitem) + "_" + str(alpha) + "_" + time.strftime("%Y-%m-%d-%H:%M:%S", time.localtime())+ ".cbp"
                f = open(data_path + "/"  + file_name,"x")
                f.write(descriptor +  '\n')
                f.write(str(Dalpha) +  '\n')
                for mu  in lst_mu:
                    f.write(str(mu) +  ' ')
                f.write('\n')
                for b  in lst_b:
                    f.write(str(b) +  ' ')
                f.write('\n')
                f.close()
