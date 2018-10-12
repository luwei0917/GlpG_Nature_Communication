#!/usr/bin/env python3
import os
import sys
import random
import time
from random import seed, randint
import argparse
import platform
from datetime import datetime
import imp
import glob
from time import sleep
import fileinput
import numpy as np
import pandas as pd
from small_script.variable_test import variable_test
from small_script.variable_test2 import variable_test2
import subprocess
from small_script.myFunctions import compute_theta_for_each_helix
from small_script.myFunctions import *

# Useful codes
# os.system("awk '{print $NF}' all_wham.dat > e_total")
# tr " " "\n"
# sed 1d
# sort -u -k 3
# sed -e 's/+T//'
# import re
# numbers = re.compile(r'(\d+)')
# def numericalSort(value):
#     parts = numbers.split(value)
#     parts[1::2] = map(int, parts[1::2])
#     return parts
# mypath = os.environ["PATH"]
# os.environ["PATH"] = "/home/wl45/python/bin:/home/wl45/opt:" + mypath
# my_env = os.environ.copy()

parser = argparse.ArgumentParser(description="This is my playground for current project")
parser.add_argument("-r", "--run", help="test mode",
                    action="store_true")
parser.add_argument("-s", "--see", help="test mode",
                    action="store_true")
# parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("-m", "--mode", type=int, default=0)
parser.add_argument("-d", "--day", type=str, default="someday")
parser.add_argument("-l", "--label", type=str, default="label")
parser.add_argument("-t", "--test", action="store_true", default=False)
args = parser.parse_args()

if args.test:
    do = print
else:
    do = os.system
cd = os.chdir

base_run_slurm = '''\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST

module load GCC/4.9.3 OpenMPI/1.8.8
srun /home/wl45/build/awsem_lipid_fluctuations/src/lmp_mpi -p 12x1 -in 2xov_{}.in\n'''

base_slurm = '''\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=01:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun {}\n'''

def replace(TARGET, FROM, TO):
    do("sed -i.bak 's/{}/{}/g' {}".format(FROM,TO,TARGET))



def continueRunConvertion(n=12, rerun=0, name="2xov", convert_read_data=False):
    rerun_plus_one = rerun + 1
    do(f"cp {name}_0.in {name}_{rerun_plus_one}.in")
    fileName = f"{name}_{rerun_plus_one}.in"
    replace(fileName, "variable r world "+ " ".join(str(i) for i in list(range(n))), "")
    replace(fileName, "# read_restart restart.25000000", "variable r world "+ " ".join(str(i) for i in list(range(n))))
    initial_steps = 20000000 * rerun_plus_one
    replace(fileName, "read_restart restart.extended", f"read_restart restart.$r.{initial_steps}")
    replace(fileName, "read_restart restart.native_topology", f"read_restart restart.$r.{initial_steps}")
    if convert_read_data:
        replace(fileName, "read_data data.2xov", f"read_restart restart.$r.{initial_steps}")
    replace(fileName, "0\/", f"{rerun_plus_one}\/")
    cmd = 'tail -n 1 log.lammps | cut -d" " -f2-'
    line = getFromTerminal(cmd).rstrip()
    replace(fileName, "reset_timestep	0", "variable w world " + line)
    replace(fileName, "fix     xbias all colvars colvars.x output x", "fix     xbias all colvars colvars.x output x.$r")
    cmd = f'grep "temper" {name}_{rerun_plus_one}.in'
    line = getFromTerminal(cmd).rstrip()
    replace(fileName, line, line + " $w")




def rerun(extra="", mode=0, useNewTable=False):
    print(mode)
    for i in range(12):
        name = f"{extra}_{i}"
        do(f"mkdir -p {name}")
        cd(str(name))
        do("cp ~/opt/2xov_eval/* .")
        if mode == 0:
            do("cp 2xov_eval_mar31.in 2xov_eval.in")
            do(f"cp 2xov_{extra}.seq 2xov.seq")
        elif mode == 1:
            do("cp 2xov_eval_apr13.in 2xov_eval.in")
        if mode == 2:
            do("cp 2xov_eval_n_side.in 2xov_eval.in")
        if mode == 3:
            do("cp 2xov_eval_n_side_3_helix.in 2xov_eval.in")
        if mode == 4:
            do("cp 2xov_eval_n_side_4_helix.in 2xov_eval.in")
        if useNewTable:
            create_zim("2xov.seq", isNew=False)
        else:
            create_zim("2xov.seq", isNew=True)
        fileName = "2xov_eval.in"
        if mode == 1:  # change lipid environment
            replace(fileName, "DMPC", extra)
        replace(fileName, "MY_RUN", str(i))
        # do("cp ../../zim .")
        replace("run.slurm", "ctbp-common", "commons")
        # replace("run.slurm", "ctbp-common", "interactive")
        do("sbatch run.slurm")
        cd("..")

# def extra(fileName, offset=0):
#     replace(fileName, "1 1 30 5", "1 1 30 {}".format(offset))
# def rerun(extra=extra, offset=0):
#     do("cp 2xov_0.in rerun_{}.in".format(offset))
#     fileName = "rerun_{}.in".format(offset)
#     replace(fileName, "fix               1 all nve", "")
#     replace(fileName, "fix               2 all langevin", "#")
#     replace(fileName, "0\/", "recompute_offset_{}\/".format(offset))
#     replace(fileName, "dump		1 all atom 4000", "#")
#     replace(fileName, "dump_modify	1 sort id", "")
#     replace(fileName, "run", "#")
#     replace(fileName, "restart         100000 restart", "rerun 0\/dump.lammpstrj dump x y z")

#     # replace(fileName, "1 1 30 5", "1 1 30 0")
#     # extra(fileName, offset)
#     slurm = "rerun_{}.slurm".format(offset)
#     do("cp ~/opt/2xov_eval/run.slurm " + slurm)
#     replace(slurm, "2xov_eval", "rerun_{}".format(offset))
#     # replace(slurm, "commons", "ctbp-common")
#     do("mkdir recompute_offset_{}".format(offset))
#     do("sbatch " + slurm)

def scancel_jobs_in_folder(folder):
    cd(bias)
    cmd = "find -name 'slurm-*' | rev | awk -F'[-.]' '{print $2}' | rev"
    lines = getFromTerminal(cmd).splitlines()
    for line in lines:
        # print(line)
        do("scancel " + line)
    cd("..")

localQ_slurm = '''#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:30:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun python3 ~/opt/gg_server.py -d feb28 -m 2
'''
quick_slurm = '''#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:30:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun python3 ~/opt/gg_server.py -d mar03 -m 3
'''

freeEnergy = """\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=30G
#SBATCH --time=23:00:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun python2 ~/opt/pulling_compute-pmf.py {}
"""
quick_template_slurm = '''\
#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=01:30:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun {}
'''
quick_template_large_mem_slurm = '''#!/bin/bash
#SBATCH --job-name=CTBP_WL
#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=01:30:00
#SBATCH --mail-user=luwei0917@gmail.com
#SBATCH --mail-type=FAIL
echo "My job ran on:"
echo $SLURM_NODELIST
srun python3 ~/opt/gg_server.py {}
'''

def compute_quantity(cmd="", queue=1, sim_list=["0"], bias=""):
    """Compute some quantity.
    Queue 0 is ctbp-common, 1 is interactive, 2 is commons.
    bias is pre name of simulation folder
    """
    # print(cmd)
    simulation_list = glob.glob(f"{bias}*")
    print(simulation_list)
    dateAndTime = datetime.today().strftime('%d_%h_%H%M%S')
    # print(sim_list)
    for sim in sim_list:
        for folder in simulation_list:
            cd(folder)
            cd(sim)
            print(folder)
            quick = quick_template_slurm.format(cmd)
            with open(f"quick_{dateAndTime}.slurm", "w") as f:
                if queue == 1:
                    quick = quick.replace("--time=01:30:00", "--time=00:30:00")
                    quick = quick.replace("#SBATCH --account=ctbp-common", "")
                    quick = quick.replace("ctbp-common", "interactive")
                if queue == 2:
                    quick = quick.replace("ctbp-common", "commons")
                f.write(quick)
            do(f"sbatch quick_{dateAndTime}.slurm")
            cd("../..")

def compute_completeZ(temper=False, **kwargs):
    print("compute completeZ")
    if temper:
        cmd = "python3 ~/opt/gg_server.py -d mar10 -m 7"
    else:
        cmd = "python3 ~/opt/gg_server.py -d mar16 -m 3"
    compute_quantity(cmd=cmd, **kwargs)

def compute_disReal(temper=False, targetMode=0, name="2xov", **kwargs):
    print("compute DisReal")
    # if targetMode == 0:
    if temper:
        cmd = f"python3 ~/opt/small_script/temper_compute_HQ.py {name} -t {targetMode}"
    else:
        cmd = f"python3 ~/opt/small_script/find_distance.py {name} -t {targetMode}"
    # elif targetMode == 1:
    #     if temper:
    #         cmd = "python3 ~/opt/gg_server.py -d apr01 -m 3"
    #     else:
    #         cmd = f"python3 ~/opt/small_script/find_distance.py -t {targetMode}"
    compute_quantity(cmd=cmd, **kwargs)

# def compute_NsideEnergy(temper=False, targetMode=0, name="2xov", **kwargs):
#     print("compute DisReal")
#     # if targetMode == 0:
#     if temper:
#         cmd = f"python3 ~/opt/gg_server.py -d mar17 -m 1"
#     else:
#         cmd = f"python3 ~/opt/small_script/find_distance.py {name} -t {targetMode}"
#     # elif targetMode == 1:
#     #     if temper:
#     #         cmd = "python3 ~/opt/gg_server.py -d apr01 -m 3"
#     #     else:
#     #         cmd = f"python3 ~/opt/small_script/find_distance.py -t {targetMode}"
#     compute_quantity(cmd=cmd, **kwargs)

def let_compute_localQ(temper=False, **kwargs):
    print("compute LocalQ")
    if temper:
        cmd = "python3 ~/opt/gg_server.py -d mar28 -m 1"
    # if temper:
    #     cmd = "python3 ~/opt/gg_server.py -d feb28 -m 2"
    # else:
    #     pass
        # cmd = "python3 ~/opt/small_script/find_distance.py"
        #         native_contacts_table = compute_localQ_init()
        # for i in range(12):
        #     compute_localQ(native_contacts_table, pre=".", ii=i)
    compute_quantity(cmd=cmd, **kwargs)

if args.day == "common":
    if args.mode == 4:
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # simulation_list = ['dis_30.0']
        # simulation_list = ['dis_116.0', 'dis_332.0', 'dis_128.0', 'dis_266.0', 'dis_296.0', 'dis_290.0', 'dis_314.0', 'dis_176.0', 'dis_272.0', 'dis_284.0', 'dis_158.0', 'dis_338.0']
        # simulation_list = ['dis_326.0', 'dis_206.0', 'dis_254.0', 'dis_344.0', 'dis_308.0', 'dis_134.0', 'dis_152.0', 'dis_194.0', 'dis_320.0', 'dis_200.0', 'dis_212.0', 'dis_110.0', 'dis_248.0', 'dis_188.0', 'dis_242.0', 'dis_218.0', 'dis_350.0', 'dis_164.0', 'dis_236.0', 'dis_146.0', 'dis_182.0', 'dis_140.0', 'dis_122.0', 'dis_302.0', 'dis_224.0', 'dis_230.0', 'dis_278.0', 'dis_260.0', 'dis_170.0']
        # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            # cd("0")
            cd("1")
            # cd("2")
            # cd("3")
            # rerun(extra="Go", mode=2)
            # rerun(extra="Go_3helix", mode=3)
            rerun(extra="Go_4helix", mode=4)
            cd("../..")
    if args.mode == 3:
        goEnergy = False
        goEnergy3H = False
        goEnergy4H = True
        rerun = 3
        end = 2
        cwd = os.getcwd()
        print(cwd)
        pre = '/'.join(cwd.split("/")[:-2]) + "/"
        print(pre)
        # exit()
        # pre = "/scratch/wl45/apr_2018/sixth/"
        data_folder = "/scratch/wl45/aug_2018/02_week/freeEnergy/all_data_folder/"
        # folder_list = ["rg_0.15_lipid_1.0_mem_1_go_0.8_long"]
        folder_list = [cwd.split("/")[-2]]
        # label = "sixth_long"
        label = args.label
        # with open(label, "w") as f:
        #     f.write("test\n")
        # cd("simulation")
        # exit()
        # folder_list = ["23oct/memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rgWidth_memb_3_rg_0.1_lipid_1_extended",
        #                 "rgWidth_memb_3_rg_0.1_lipid_1_topology",
        #                 "expand_distance_rgWidth_memb_3_rg_0.1_lipid_1_extended"]

        for attempt in range(2):
            try:
                process_complete_temper_data_3(pre, data_folder, folder_list, rerun=rerun, end=end, average_z=True, disReal=True, dis_h56=True, localQ=False, goEnergy=goEnergy, goEnergy3H=goEnergy3H, goEnergy4H=goEnergy4H, label=label)
                break
            except FileNotFoundError:
                bias = "dis"
                simulation_list = glob.glob(f"{bias}_*")
                # simulation_list = ['dis_30.0']
                # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
                print(simulation_list)
                for dis in simulation_list:
                    print(dis)
                    cd(dis)
                    i = rerun
                    i_plus_one = i +1
                    if os.path.exists(f"log{i}"):
                        do(f"mv log{i} back_log{i}")  # in case of override
                    do(f"mkdir -p log{i}")
                    do(f"cp log.lammps log{i}/")
                    cd("..")
    if args.mode == 2:
        bias = "dis"
        simulation_list = glob.glob(f"{bias}_*")
        # simulation_list = ['dis_30.0']
        # simulation_list = ['dis_86.0', 'dis_84.0', 'dis_76.0', 'dis_72.0', 'dis_54.0', 'dis_70.0', 'dis_50.0', 'dis_56.0', 'dis_80.0', 'dis_30.0', 'dis_88.0', 'dis_44.0', 'dis_46.0', 'dis_96.0', 'dis_38.0']
        print(simulation_list)
        for dis in simulation_list:
            print(dis)
            cd(dis)
            i = 2
            i_plus_one = i +1
            # do(f"mv log{i} back_log{i}")  # in case of override
            # do(f"mkdir -p log{i}")
            # do(f"cp log.lammps log{i}/")

            # my_file = sys.Path(f"{i_plus_one}")
            # if my_file.is_dir():
            #     print("Attension")
            #     exit()
            continueRunConvertion(n=12, rerun=i)
            # continueRunConvertion(n=12, rerun=i, convert_read_data=True)
            do(f"mkdir {i_plus_one}")
            # do(f"sed 's/2xov_{i}/2xov_{i_plus_one}/g' run_{i}.slurm > run_{i_plus_one}.slurm")
            # replace(f"run_{i_plus_one}.slurm", "/home/ns24/lmp_mpi", "/home/wl45/build/awsem_lipid_fluctuations/src/lmp_mpi")
            # do(f"sbatch run_{i_plus_one}.slurm")
            run_slurm = base_run_slurm.format(i_plus_one)
            with open(f"run_{i_plus_one}.slurm", "w") as r:
                r.write(run_slurm)
            # replace(f"run_{i_plus_one}.slurm", "ctbp-common", "commons")
            do(f"sbatch run_{i_plus_one}.slurm")
            cd("..")
    if args.mode == 1:
        queue = 2
        # i = "0"
        # i = "1"
        i = "2"
        # i = "3"
        # let_compute_localQ(temper=True, bias="dis_", sim_list=[i], queue=1)

        # compute_disReal(temper=True, bias="dis_", sim_list=[i], queue=queue)
        # compute_disReal(temper=True, targetMode=1, bias="dis_", sim_list=[i], queue=queue)
        # compute_completeZ(temper=True, bias="dis_", sim_list=[i], queue=queue)

        compute_disReal(temper=True, targetMode=2, bias="dis_", sim_list=[i], queue=queue)
        compute_disReal(temper=True, targetMode=3, bias="dis_", sim_list=[i], queue=queue)

if args.day == "freeEnergyCalculation":
    if args.mode == 1:
        temp_list = ["all"]
        data_folder = "all_data_folder/"
        bias_list = {"2d_z_qw":"13", "1d_dis":"9", "2d_z_dis":"14", "2d_qw_dis":"11", "1d_qw":"10", "1d_z":"12"}
        bias_list = {"2d_zAverage_dis":"17", "2d_z_qw":"13"}

        # bias_list = {"2d_z_qw_enhance":"16", "2d_z_qw":"13", "1d_dis":"9"}  # z and Dis_h56
        i = -2
        # freeEnergy_folder = f"sixth_i235d/"
        # freeEnergy_folder = f"sixth_i255d/"
        # freeEnergy_folder = f"sixth_long/"
        freeEnergy_folder = f"second_start_extended_combined"
        print(freeEnergy_folder)
        # folder_list = ["memb_3_rg_0.1_lipid_1_extended"]
        # folder_list = ["rerun_1_08_Mar_154259"]
        # folder_list = [f"first_rerun_{sample_range_mode}_12_Mar_151630" for i in range(4,6)]
        # folder_list = [f"sixth_i235drerun_3_03_Apr_220358"]
        # folder_list = [f"sixth_i255drerun_3_04_Apr_145735"]
        folder_list = [f"second_start_extended_combined_may19"]
        # submode_list = ["_no_energy"]
        # submode_list = ["", "only_500"]
        # submode_list = ["350", "400", "450", "500", "550"]

        # temp_dic = {"_350-550":["350", "400", "450", "500", "550"]}
        # temp_dic = {"_280-350":["300", "335", "373"]}
        temp_dic = {"_280-350":["335", "373", "417"]}
        temp_dic = {"_280-350":["335", "373"]}
        # temp_dic = {"_280-350":["280", "290", "300", "315", "335"]}
        # dic = {"T0":280, "T1":300, "T2":320, "T3":350, "T4":375, "T5":400
        # temp_dic = {"_280-350":["280", "290", "300", "310", "320", "335", "350"]}
        for temp_mode, temp_list in temp_dic.items():
            move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=4, sample_range_mode=i, sub_mode_name=temp_mode, average_z=5, chosen_mode=9)  # chosen_mode 4 use Dis_h56

        cd(freeEnergy_folder)
        for temp_mode, temp_list in temp_dic.items():
                cd(temp_mode)
                for bias, mode in bias_list.items():
                    # name = "low_t_" + bias
                    name = str(bias)
                    print(name)
                    do("rm -r "+name)
                    do("mkdir -p " + name)
                    cd(name)
                    make_metadata_3(temps_list=temp_list,k=0.02, i=i)
                    # nsample = len(folder_list)*2500
                    nsample = len(folder_list)*5000
                    do(f"python3 ~/opt/pulling_analysis_2.py -m {mode} --commons 1 --nsample {nsample} --submode 21 --force 8")
                    cd("..")
                cd("..")
        cd("..")


def pick_structure_generate_show_script(n=2):
    with open("show.pml", "w") as f:
        for structure_index in range(0, n):
            f.write("load structure_%s.pdb\n" % structure_index)
            f.write("cealign structure_0, structure_%s\n" % structure_index)
            f.write("spectrum count, rainbow_rev, structure_%s, byres=1\n" % structure_index)
        # f.write("hide lines, all\n")
        # f.write("show cartoon, all\n")
        # f.write("hide nonbonded, all\n")

if args.day == "jun23":
    if args.mode == 2:
        cmd_pre = "python2 ~/opt/script/BuildAllAtomsFromLammps.py"
        # location_pre = "/Users/weilu/Research/server/apr_2018/sixth/rg_0.15_lipid_1.0_mem_1_go_0.8/simulation"

        # location_pre = "/Users/weilu/Research/server/may_2018/second_long/simulation"
        # cmd = cmd_pre + " " + location + " structure_2 4080 -seq ~/opt/pulling/2xov.seq"
        pick_list = ["low_e_jun23"]
        for picked in pick_list:
            do(f"mkdir {picked}")
            cd(picked)
            tt = pd.read_csv(f"/scratch/wl45/jun_2018/{picked}.csv", index_col=0)
            sample = tt.reset_index(drop=True)
            # sample["Frame"] = ((sample["Step"] - 2e7*rerun)/4000).astype("int")
            sample["rerun"] = (sample["Step"] // 2e7).astype(int)
            sample["Frame"] = ((sample["Step"] % 2e7)/4000).astype("int")
            for index, row in sample.iterrows():
                BiasTo = row["BiasTo"]
                Run = row["Run"]
                Frame = row["Frame"]
                rerun = row["rerun"]
                print(BiasTo, Run, Frame)

                # try:
                location_pre = "/scratch/wl45/may_2018/second/simulation"
                location = location_pre + f"/dis_{BiasTo}/{rerun}/dump.lammpstrj.{int(Run)}"
                cmd = cmd_pre + " " + location + f" structure_{index} {int(Frame)} -seq ~/opt/pulling/2xov.seq"
                hasProblem = do(cmd)
                if hasProblem == 0:
                    continue
                print(do(cmd))
                    # print(cmd)
                # except IOError:
                # print("-----------hi------------")
                location_pre = "/scratch/wl45/may_2018/second_long/simulation"
                location = location_pre + f"/dis_{BiasTo}/{rerun}/dump.lammpstrj.{int(Run)}"
                cmd2 = cmd_pre + " " + location + f" structure_{index} {int(Frame)} -seq ~/opt/pulling/2xov.seq"
                do(cmd2)
                    # print(cmd2)
            pick_structure_generate_show_script(n=len(sample))
            cd("..")

if args.day == "jun17":
    if args.mode == 2:
        cmd_pre = "python2 ~/opt/script/BuildAllAtomsFromLammps.py"
        # location_pre = "/Users/weilu/Research/server/apr_2018/sixth/rg_0.15_lipid_1.0_mem_1_go_0.8/simulation"

        # location_pre = "/Users/weilu/Research/server/may_2018/second_long/simulation"
        # cmd = cmd_pre + " " + location + " structure_2 4080 -seq ~/opt/pulling/2xov.seq"
        pick_list = ["low_e_jun01_h56", "low_e_jun01_h34", "low_e_jun01_h12", "low_e_jun01_out", "low_e_jun01_pre",
                        "low_e_jun01_transition", "low_e_jun01_post_transition"]
        pick_list = ["low_e_path1", "low_e_path2"]
        for picked in pick_list:
            do(f"mkdir {picked}")
            cd(picked)
            tt = pd.read_csv(f"/scratch/wl45/jun_2018/{picked}.csv", index_col=0)
            sample = tt.reset_index(drop=True)
            # sample["Frame"] = ((sample["Step"] - 2e7*rerun)/4000).astype("int")
            sample["rerun"] = (sample["Step"] // 2e7).astype(int)
            sample["Frame"] = ((sample["Step"] % 2e7)/4000).astype("int")
            for index, row in sample.iterrows():
                BiasTo = row["BiasTo"]
                Run = row["Run"]
                Frame = row["Frame"]
                rerun = row["rerun"]
                print(BiasTo, Run, Frame)

                # try:
                location_pre = "/scratch/wl45/may_2018/second/simulation"
                location = location_pre + f"/dis_{BiasTo}/{rerun}/dump.lammpstrj.{int(Run)}"
                cmd = cmd_pre + " " + location + f" structure_{index} {int(Frame)} -seq ~/opt/pulling/2xov.seq"
                hasProblem = do(cmd)
                if hasProblem == 0:
                    continue
                print(do(cmd))
                    # print(cmd)
                # except IOError:
                # print("-----------hi------------")
                location_pre = "/scratch/wl45/may_2018/second_long/simulation"
                location = location_pre + f"/dis_{BiasTo}/{rerun}/dump.lammpstrj.{int(Run)}"
                cmd2 = cmd_pre + " " + location + f" structure_{index} {int(Frame)} -seq ~/opt/pulling/2xov.seq"
                do(cmd2)
                    # print(cmd2)
            pick_structure_generate_show_script(n=len(sample))
            cd("..")

        # tt = pd.read_csv("/scratch/wl45/jun_2018/low_e_jun01_h56.csv", index_col=0)
        # tt = pd.read_csv("/scratch/wl45/jun_2018/low_e_jun01_h34.csv", index_col=0)
        # tt = pd.read_csv("/scratch/wl45/jun_2018/low_e_jun01_h12.csv", index_col=0)
        # tt = pd.read_csv("/scratch/wl45/jun_2018/low_e_jun01_out.csv", index_col=0)
        # tt = pd.read_csv("/scratch/wl45/jun_2018/low_e_jun01_pre.csv", index_col=0)
        # tt = pd.read_csv("/scratch/wl45/jun_2018/low_e_jun01_transition.csv", index_col=0)
        # tt = pd.read_csv("/scratch/wl45/jun_2018/low_e_jun01_post_transition.csv", index_col=0)
        # # rerun = 1
        # # sample = tt.sample(5).reset_index(drop=True)
        # ample = tt.reset_index(drop=True)
        # # sample["Frame"] = ((sample["Step"] - 2e7*rerun)/4000).astype("int")
        # sample["rerun"] = (sample["Step"] // 2e7).astype(int)
        # sample["Frame"] = ((sample["Step"] % 2e7)/4000).astype("int")
        # for index, row in sample.iterrows():
        #     BiasTo = row["BiasTo"]
        #     Run = row["Run"]
        #     Frame = row["Frame"]
        #     rerun = row["rerun"]
        #     print(BiasTo, Run, Frame)

        #     location = location_pre + f"/dis_{BiasTo}/{rerun}/dump.lammpstrj.{int(Run)}"
        #     cmd = cmd_pre + " " + location + f" structure_{index} {int(Frame)} -seq ~/opt/pulling/2xov.seq"
        #     print(cmd)
        #     do(cmd)
        # pick_structure_generate_show_script(n=len(sample))
    if args.mode == 1:
        data = pd.read_feather("/scratch/wl45/may_2018/03_week/all_data_folder/second_start_extended_combined_may19.feather")
        data = data.reset_index(drop=True)
        # data["BiasedEnergy"] = data["TotalE"] + 0.2*data["AMH_4H"]
        data["BiasedEnergy"] = data["Lipid"] + data["Rg"] + data["Membrane"] + data["AMH-Go"] + 0.2*data["AMH_4H"]
        data["BiasEnergy"] = 0.02 * (data["BiasTo"] - data["DisReal"])**2
        data["Energy_with_all_bias"] = data["BiasEnergy"] + data["BiasedEnergy"]

        t_pos = data.query("TempT == 373 and DisReal > 52 and DisReal < 57 and z_average > -4 and z_average < 0").reset_index(drop=True)
        chosen = t_pos.query("Lipid1 < -0.5").sort_values("Energy_with_all_bias").head(n=20)
        # chosen = t_pos.sort_values("Energy_with_all_bias").head(n=20)
        chosen.to_csv("low_e_jun01_pre.csv")

        t_pos = data.query("TempT == 373 and DisReal > 57 and DisReal <63 and z_average > -5 and z_average < -2").reset_index(drop=True)
        chosen = t_pos.query("Lipid1 < -0.5 and Lipid10 < -0.5").sort_values("Energy_with_all_bias").head(n=20)
        # chosen = t_pos.sort_values("Energy_with_all_bias").head(n=20)
        chosen.to_csv("low_e_jun01_transition.csv")

        t_pos = data.query("TempT == 373 and DisReal > 63 and DisReal <72 and z_average > -6 and z_average < -3").reset_index(drop=True)
        chosen = t_pos.query("Lipid1 < -0.5").sort_values("Energy_with_all_bias").head(n=20)
        # chosen = t_pos.sort_values("Energy_with_all_bias").head(n=20)
        chosen.to_csv("low_e_jun01_post_transition.csv")

        t_pos = data.query("TempT == 373 and DisReal > 80 and DisReal < 100 and z_average > -8 and z_average < -4").reset_index(drop=True)
        chosen = t_pos.query("Lipid1 < -0.5").sort_values("Energy_with_all_bias").head(n=20)
        chosen.to_csv("low_e_jun01_h56.csv")

        t_pos = data.query("TempT == 373 and DisReal > 140 and DisReal < 180 and z_average > -14 and z_average < -8").reset_index(drop=True)
        chosen = t_pos.sort_values("Energy_with_all_bias").head(n=20)
        chosen.to_csv("low_e_jun01_h34.csv")

        t_pos = data.query("TempT == 373 and DisReal > 220 and DisReal < 250 and z_average > -14 and z_average < -10").reset_index(drop=True)
        chosen = t_pos.sort_values("Energy_with_all_bias").head(n=20)
        chosen.to_csv("low_e_jun01_h12.csv")

        t_pos = data.query("TempT == 373 and DisReal > 260 and z_average < -16").reset_index(drop=True)
        chosen = t_pos.sort_values("Energy_with_all_bias").head(n=20)
        chosen.to_csv("low_e_jun01_out.csv")
