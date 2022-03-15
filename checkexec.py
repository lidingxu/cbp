from pathlib import Path
import os
import sys

# 1: solved. 2: unsolved
def main(argv):
    test_path = os.getcwd() + "/test"
    path = os.getcwd()
    benckmarkname=argv[1] 
    solvername=argv[2]
    dataname=argv[3]
    log_path = test_path + "/" + benckmarkname + "/logs/" + solvername + "/" + dataname +  "_" + solvername + ".log"
    result_path = test_path + "/" + benckmarkname + "/results/"+ solvername + "/" + dataname +  "_" + solvername + ".result"
    #print(benckmarkname, solvername, dataname, log_path, result_path)
    if solvername == "bcsocp":
        if os.path.isfile(result_path) and os.path.getsize(result_path) > 0:
            return 1
        else:
            if os.path.isfile(log_path):
                os.remove(log_path)
            if os.path.isfile(result_path):
                os.remove(result_path)
            return 2
    else:
        if os.path.isfile(log_path) and os.path.getsize(log_path) > 0:
            return 1
        else:
            if os.path.isfile(log_path):
                os.remove(log_path)
            return 2
    return 0

if __name__ == "__main__":
    rtcode = main(sys.argv)
    sys.exit(rtcode)