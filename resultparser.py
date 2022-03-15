import os
import csv
from typing import NamedTuple
import collections
import numpy as np
from pathlib import Path
import os

max_primal_bound = 2**31
min_dual_bound = -1
time_limit = 3600
solvers =  [ "bphybridph", "bphybridphad3"]
benchmarks = ["CloudMedium"]




def extract_scip_csv(file_path, file):
    ls = open(file_path).readlines()
    #print(file)
    stat_keys = ["Total Time", "CKNAP_Pricer", "nodes (total)", "Primal Bound", "First Solution", "Dual Bound", "Gap", "pricing column exact", "pricing log sum shifted gap"]
    stat_dict = {}
    entries = {}
    for l in ls:
        for stat_key in stat_keys:
            if  l.split(":")[0].strip() == stat_key:
                stat_dict[stat_key] = l
    entries["total_time"] = float(stat_dict["Total Time"].split()[3])
    if entries["total_time"] > time_limit * 1.5:
        print(file_path)
    entries["#nodes"] = int(stat_dict["nodes (total)"].split()[3])
    entries["primal_bound"] = float(stat_dict["Primal Bound"].split()[3])
    entries["dual_bound"] = float(stat_dict["Dual Bound"].split()[3])
    entries["relative_gap(%)"] = abs(entries["dual_bound"] - entries["primal_bound"])/ entries["primal_bound"] * 100 #float(stat_dict["Gap"].split()[2])
    entries["pricing_time"]= float(stat_dict["CKNAP_Pricer"].split()[2])
    entries["pricing_time(%)"]= float(stat_dict["CKNAP_Pricer"].split()[2]) / float(stat_dict["Total Time"].split()[3]) * 100
    entries["#pricing_calls"] =  int(stat_dict["CKNAP_Pricer"].split()[4])
    entries["#columns_gen"] =  int(stat_dict["CKNAP_Pricer"].split()[5])
    entries["columns_gen_exact"] =  float(stat_dict["pricing column exact"].split(":")[1]) / entries["#columns_gen"] * 100
    entries["#columns_gen_exact"] =  float(stat_dict["pricing column exact"].split(":")[1]) 
    entries["logsum shift pricing gap"] =  float(stat_dict["pricing log sum shifted gap"].split(":")[1])
    entries["absolute_gap"] =  float(entries["primal_bound"] ) - entries["dual_bound"] 
    entries["first_solution"] = float(stat_dict["First Solution"].split()[3]) 
    entries["primal_bound close"] = abs(float(stat_dict["First Solution"].split()[3]) - entries["primal_bound"]) /  float(stat_dict["First Solution"].split()[3]) * 100
    entries["primal_bound improved"] =  abs(entries["first_solution"] -  entries["primal_bound"])
    return entries

def extract_cplex_csv(file_path, file):
    ls = open(file_path).readlines()
    #print(file)
    stat_keys = ["soltime", "nodes", "relgap", "obj", "best_bound", "first_solution"]
    stat_dict = {}
    entries = {}
    for l in ls:
        for stat_key in stat_keys:
            if  l.split(":")[0].strip() == stat_key:
                stat_dict[stat_key] = l
    #print(stat_dict["soltime"].split(),  stat_dict["nodes"].split(), stat_dict["obj"].split(), stat_dict["relgap"].split())
    entries["total_time"] = float(stat_dict["soltime"].split()[1])
    if entries["total_time"] > time_limit * 1.5:
        print(file_path)
    entries["#nodes"] = int(stat_dict["nodes"].split()[1])
    entries["primal_bound"] = float(stat_dict["obj"].split()[1])
    entries["dual_bound"] = float(stat_dict["best_bound"].split()[1])
    entries["first_solution"] = entries["primal_bound"]
    #entries["first_solution"] = float(stat_dict["first_solution"].split()[1])
    if( -entries["dual_bound"] + entries["primal_bound"] < 1 + 1e-6) :
         entries["dual_bound"] = entries["primal_bound"]
    entries["relative_gap(%)"] = abs(entries["dual_bound"] - entries["primal_bound"])/ entries["primal_bound"] * 100# float(stat_dict["relgap"].split()[1])*10000
    entries["absolute_gap"] =  float(entries["primal_bound"] ) - entries["dual_bound"]  
    entries["primal_bound close"] = abs(entries["primal_bound"] - entries["first_solution"]) /  entries["first_solution"] * 100
    entries["primal_bound improved"] =  abs(entries["first_solution"] -  entries["primal_bound"])
    return entries
    #writer.writerow([name, algo, primal_bound, gap, total_time, pricing_time, pricing_calls, columns_gen, nodes])

def extract_sol(file_path, file):
    ls = open(file_path).readlines()
    stat_keys = [ "primal_bound", "dual_bound"]
    stat_dict = {}
    entries = {}
    for l in ls:
        for stat_key in stat_keys:
            if  l.split(":")[0].strip() == stat_key:
                stat_dict[stat_key] = l
    entries["primal_bound"] = float(stat_dict["primal_bound"].split()[1])
    entries["dual_bound"] = float(stat_dict["dual_bound"].split()[1])
    return entries

def get_method(mystr):
    if mystr == "g":
        return "gaussian case"
    elif mystr == "h":
        return "Hoeffding's Inequality"
    else:
        return "Distributionally robust formulation"

def extract_data_csv(file_path, file):
    ls = open(file_path).readlines()
    entries = {}
    descriptor = ls[0].split(" ")
    entries["uncertainty type"] = get_method(descriptor[0])
    entries["capacity"] = float(descriptor[1])
    entries["#items"]= float(descriptor[2])
    entries["confidence level"]= float(descriptor[3])
    entries["Dalpha"]= float(ls[1])
    avga =  0
    avgb =  0
    numitem = int(descriptor[2])
    for a in ls[2].split(" "):
        if a != "\n":
            avga += float(a)
    for b in ls[3].split(" "):
        if b != "\n":
            avgb += float(b)  
    entries["avg. a"] = round(avga / numitem, 3)
    entries["avg. b"] = round(avgb / numitem, 3)
    return entries


shift = { "relative_gap(%)": 1.0,"total_time": 1.0, "total": 0, "#nodes": 1.0, "#columns_gen": 1.0, "columns_gen_exact": 1.0, "logsum shift pricing gap":1.0, "primal_bound close":1.0, "primal_bound improved":1.0, "pricing_time(%)": 1.0}

def Stat(name):
    return {"algorithm": name, "solved": 0, "improved": 0,  "relative_gap(%)": 0.0, "total_time":0.0, "#nodes":  0.0,  "#columns_gen":0.0, "columns_gen_exact":0.0, "#columns_gen_exact": 0.0, "logsum shift pricing gap":0.0, "primal_bound close":0.0, "primal_bound improved":0.0, "total": 0, "pricing_time(%)": 0} 


def add(stat, entries):
    stat["solved"] += 1 if  float(entries["absolute_gap"]) < 0.99 else 0
    stat["improved"] += 1 if  float(entries["primal_bound"])  < float(entries["first_solution"]) - 0.9 else 0
    stat["total"] += 1
    keys = ["relative_gap(%)", "total_time", "#nodes", "columns_gen_exact", "#columns_gen",  "primal_bound close", "pricing_time(%)", "primal_bound improved"]
    for key in keys:
        if key in entries:
            stat[key] += np.log(float(entries[key])+ shift[key])
    keys = ["logsum shift pricing gap", "#columns_gen_exact"]
    for key in keys:
        if key in entries:
            stat[key] += float(entries[key])


def avgStat(stat, solver):
    display_keys = [ "total_time", "relative_gap(%)", "primal_bound close", "#nodes", "solved", "improved"]
    if solver  == "bcsocp":
        keys = ["relative_gap(%)", "total_time", "#nodes",   "primal_bound close",  "primal_bound improved"]
        for key in keys:
            stat[key] = np.exp(stat[key] / stat["total"])- shift[key]
    else:
        keys = ["relative_gap(%)", "total_time", "#nodes", "#columns_gen",  "columns_gen_exact", "primal_bound close",  "primal_bound improved", "pricing_time(%)"]
        for key in keys:
            stat[key] = np.exp(float(stat[key] / stat["total"] ))- shift[key]
        key = "logsum shift pricing gap"
        stat[key] /= stat["#columns_gen_exact"]
        stat[key] = np.exp(float(stat[key])) - shift[key]
        display_keys.extend(["#columns_gen", "columns_gen_exact", "logsum shift pricing gap", "pricing_time(%)"])
    avg_info = {}
    display = ""
    for key in display_keys:
        avg_info[key] = stat[key]
        if key == "solved" or key == "improved":
            display += str(round(stat[key],3)) + "/"+  str(stat["total"]) + " & "
        else:
            display += str(round(stat[key],3)) + " & "
    print(avg_info , solver)
    print(display)


for benchmark in benchmarks:
    # parse all results and logs, create solution for each instance
    data_dir_path = os.getcwd() + "/data/" + benchmark
    solution_dir_path = os.getcwd() + "/test/"+benchmark+"/solutions"
    Path(solution_dir_path).mkdir(parents=True, exist_ok=True)
    instances = os.listdir(data_dir_path)
    for instance in instances:
        len_name = len(instance)
        entries_sol = {"primal_bound": max_primal_bound, "dual_bound": min_dual_bound}
        for solver in solvers:
            if solver == "bcsocp":
                out_dir_path = os.getcwd() + "/test/" + benchmark + "/results/" + solver 
            else:
                out_dir_path = os.getcwd() + "/test/" + benchmark + "/logs/" + solver
            out_files = os.listdir(out_dir_path)
            for out_file in out_files:
                if instance == out_file[0:len_name]:
                    if solver == "bcsocp":
                        entries = extract_cplex_csv(out_dir_path  + "/" + out_file ,  out_file)                
                    else:
                        entries = extract_scip_csv(out_dir_path  + "/" + out_file ,  out_file)
                    break
            entries_sol["primal_bound"] = min(entries_sol["primal_bound"],  entries["primal_bound"])
            entries_sol["dual_bound"] = max(entries_sol["dual_bound"],  entries["dual_bound"])
            solution_path = solution_dir_path + "/" + instance + ".sol"
            solution_file = open(solution_path, "w")
            solution_file.write('primal_bound: ' + str(entries_sol["primal_bound"] ) + '\n')
            solution_file.write('dual_bound: ' + str(entries_sol["dual_bound"]) + '\n')
            solution_file.close()
    # compute statistics
    solutions = os.listdir(solution_dir_path)
    for solver in solvers:
        if solver == "bcsocp":
            out_dir_path = os.getcwd() + "/test/" + benchmark + "/results/" + solver 
        else:
            out_dir_path = os.getcwd() + "/test/" + benchmark + "/logs/" + solver
        out_files = os.listdir(out_dir_path)
        stat = Stat(solver)
        for out_file in out_files:
            isfind = False
            for solution in solutions:
                instance_name = solution[0:-4]
                if instance_name == out_file[0: len(instance_name)]:
                    entries_sol = extract_sol(solution_dir_path + "/" + solution, solution)
                    isfind = True
                    break
            assert(isfind)
            if solver == "bcsocp":
                entries = extract_cplex_csv(out_dir_path  + "/" + out_file ,  out_file)
                entries["instance"] = out_file 
                entries["algorithm"] = solver
                entries["primal_bound close"] = (abs(entries["primal_bound"] - entries["first_solution"]) /  abs(max(entries["first_solution"] - entries_sol["dual_bound"], 1e-6))) * 100
                add(stat, entries)
                #print(entries)
            else:
                entries = extract_scip_csv(out_dir_path  + "/" + out_file ,  out_file)
                entries["instance"] = out_file 
                entries["algorithm"] = solver
                entries["primal_bound close"] = (abs(entries["primal_bound"] - entries["first_solution"]) /  abs(max(entries["first_solution"] - entries_sol["dual_bound"], 1e-6))) * 100
                add(stat, entries)
                #print(entries)
        avgStat(stat, solver)


