#!/usr/bin/env python

import sys

new_sys = []
for elem in sys.path:
    if ( "home" in elem or "work" in elem or "scratch" in elem):
        pass
    else:
        new_sys.append(elem)
sys.path = new_sys


import os
import subprocess
import time
import glob
import datetime

def TimestampMillisec64():
    return int((datetime.datetime.utcnow() - datetime.datetime(1970, 1, 1)).total_seconds() * 1000) 


tar_files = []
for elem in sys.argv[1:-1]:
    tar_files.append(elem)

dest = sys.argv[-1]



def do_tar_copy(tar_file, dest):


    tar_base = os.path.basename(tar_file)


    complete_file = os.path.join(dest, tar_base + ".complete")
    fail_file = os.path.join(dest, tar_base + ".fail")
    lock_file = os.path.join(dest, tar_base + ".lock")
    lock_glob = lock_file + "*"

    timeout = 1e3 * 5 * 60

    # Your lock must be younger than this if you wish to start
    start_timeout = 1e3 * 15



    def cmd2(command, wait=True):
        # print ""
        # print command
        the_command = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        if (not wait):
            return
        the_stuff = the_command.communicate()
        return str(the_stuff[0]),  str(the_stuff[1]), the_command.returncode


    start_time = time.time()


    while( time.time() - start_time < 3600 ):

        if (os.path.exists(complete_file)):
            print("Complete detected")
            return


        if (os.path.exists(fail_file)):
            print("Fail detected")
            sys.exit(1)


        our_pid = os.getpid()
        time_ms = TimestampMillisec64()

        locks = glob.glob(lock_glob)


        is_us = []

        for lock in sorted(locks):
            to_split = lock.replace(lock_file + "_", "")
            sp = to_split.split("_")

            this_time = int(sp[0])
            pid = int(sp[1])
            elapsed_time = time_ms - this_time

            # print(elapsed_time)
            if ( elapsed_time < timeout ):

                ok_start = elapsed_time < start_timeout

                is_us.append( pid == our_pid and ok_start)


        our_lock = lock_file + "_%i_%i"%(time_ms, our_pid) 

        if ( len(is_us) == 0 ):
            with open(our_lock, "w") as f:
                f.write("\n")

            time.sleep(5)
        else:


            if ( is_us[0] ):
                print("We're doing it!")
                out, err, code = cmd2("cd %s; tar -xzf %s"%(dest, tar_file))
                print(out)
                print(err)

                if ( code != 0 ):
                    with open(fail_file, "w") as f:
                        f.write("\n")
                else:
                    with open(complete_file, "w") as f:
                        f.write("\n")
            else:

                time.sleep(15)

    sys.exit(1)






for tar_file in tar_files:
    do_tar_copy(tar_file, dest)



