from pathlib import Path
import os

benchmarks = ["CloudSmall", "CloudMedium", "CloudBig"]


solvers = [ "bphybridphad3"]


data_path = os.getcwd() + "/data"
test_path = os.getcwd() + "/test"
benckmark_dirs = os.listdir(data_path)

clear_mode = False


#for b_dir in benckmark_dirs:
#    benchmarks.append(b_dir)

for bch in benchmarks:
    Path(os.getcwd() + "/test/"+bch).mkdir(parents=True, exist_ok=True)
    Path(os.getcwd() + "/test/"+bch+"/logs").mkdir(parents=True, exist_ok=True)
    Path(os.getcwd() + "/test/"+bch+"/results").mkdir(parents=True, exist_ok=True)
    for solver in solvers:
        result_solver_path = os.getcwd() + "/test/"+bch+"/results/"+solver
        log_solver_path = os.getcwd() + "/test/"+bch+"/logs/"+solver
        Path(result_solver_path).mkdir(parents=True, exist_ok=True)
        Path(log_solver_path).mkdir(parents=True, exist_ok=True)
        # clear 
        if clear_mode:
            for result_file in  Path(result_solver_path).iterdir():
                result_file.unlink()
            for log_file in  Path(log_solver_path).iterdir():
                log_file.unlink()
