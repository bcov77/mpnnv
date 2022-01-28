#!/usr/bin/env python

import os
import sys
import subprocess
import time
import socket



pdb_file = sys.argv[1]
number_of_iters = int(sys.argv[2])


with open(pdb_file, "rb") as f:
    pdb_data = f.read()


flags = os.O_NONBLOCK
out_read, out_write = os.pipe2(flags)
in_read, in_write = os.pipe2(flags)

f_read = in_read
f_write = out_write




def write_message(f_write, message, add_null=True):
    written = 0

    max_write = 4096
    if ( add_null ):
        message = message + bytes([0])

    while ( written < len(message) ):

        we_wrote, written = write_message_one_write(f_write, written, message)

        if ( we_wrote == 0 ):
            time.sleep(0.1)

def write_message_one_write(f_write, written, message):
    if ( written == len(message) ):
        return 0, written
    max_write = 4096

    write_upper = written + max_write
    if ( write_upper > len(message) ):
        write_upper = len(message)


    try:
        we_wrote = os.write(f_write, message[written:write_upper])
    except OSError as err:
        if err.errno == socket.errno.EAGAIN or err.errno == socket.errno.EWOULDBLOCK:
            we_wrote = 0
        else:
            raise  # something else has happened -- better reraise


    written += we_wrote

    return we_wrote, written


write_message(f_write, "0 test\n".encode("ascii"), add_null=True)

time.sleep(0.5)


command = "./mpnnv.py -pipes %s %s %s %s"%(out_read, in_write, out_write, in_read)

mpnnv = subprocess.Popen(command, shell=True, stdout=sys.stdout, stderr=sys.stderr, universal_newlines=True, close_fds=False)



os.close(out_read)
os.close(in_write)

expected_responses = {}
expected_responses["0"] = False

in_data_buf = bytes([])
out_data_buf = bytes([])

out_data_buf_written = 0


for i in range(number_of_iters):
    code = str(i + 1)

    expected_responses[code] = False

    out_data_buf += (code + " AB\n").encode("ascii")
    out_data_buf += pdb_data
    out_data_buf += "\n".encode('ascii')
    out_data_buf += bytes([0])



while True:

    last_out_written, out_data_buf_written = write_message_one_write( f_write, out_data_buf_written, out_data_buf)

    if ( last_out_written == 0 ):
        time.sleep(0.1)


    try:
        data = os.read(f_read, 1024)
    except OSError as err:
        if err.errno == socket.errno.EAGAIN or err.errno == socket.errno.EWOULDBLOCK:
            data = None
        else:
            raise  # something else has happened -- better reraise

    in_null_detected = False
    if ( not data is None ):
        if ( len(data) == 0 ):
            print("MPNNV closed pipe")
            break

        if ( len(data) > 0 ):
            if ( 0 in data ):
                in_null_detected = True
            in_data_buf += data

    if ( in_null_detected ):

        null_offset = in_data_buf.find(0)

        this_data = in_data_buf[:null_offset].decode("ascii")
        in_data_buf = in_data_buf[null_offset+1:]

        sp = this_data.strip().split()
        assert(sp[0] in expected_responses)
        assert(not expected_responses[sp[0]])
        expected_responses[sp[0]] = True

        print(this_data.strip())

    all_good = True
    for key in expected_responses:
        if ( not expected_responses[key] ):
            all_good = False
            break

    if ( all_good ):
        break


print(mpnnv.poll())

# This is what tells mpnnv to stop
os.close(f_write)
os.close(f_read)

time.sleep(0.1)

print(mpnnv.poll())

time.sleep(0.1)

print(mpnnv.poll())

print("Tester quit")







