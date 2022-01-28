#!/usr/bin/env python

import sys
import os
import argparse

sys.path = ["/home/bcov/.conda/envs/mypy_36/lib/python3.6/", "./"]

os.environ['XLA_PYTHON_CLIENT_PREALLOCATE']='false'

import generate_sequences_s2s_chain as mpnn_util
import torch
import time
import numpy as np
import socket


parser = argparse.ArgumentParser()
parser.add_argument( "-checkpoint_path", type=str, default='/projects/ml/struc2seq/data_for_complexes/training_scripts/models/new_frame_v4/checkpoints/epoch35_step175000.pt' )
parser.add_argument( "-temperature", type=float, default=0.1, help='An a3m file containing the MSA of your target' )
parser.add_argument( "-augment_eps", type=float, default=0.05, help='The variance of random noise to add to the atomic coordinates (default 0.05)' )
parser.add_argument( "-do_predictor", action="store_true", help='Whether or not to run the predictor step and generate an estimated ddg for each MPNN sequence (default False)' )
parser.add_argument( "-protein_features", type=str, default='full', help='An a3m file containing the MSA of your target' )
parser.add_argument( "-num_connections", type=int, default=30, help='Number of neighbors each residue is connected to, default 64, higher number leads to better interface design but will cost more to run the model.' )
parser.add_argument( "-pipes", type=int, nargs="*", help='read_pipe write_pipe close_pipe close_pipe' )
parser.add_argument( "-test_pdb", type=str, default="")

args = parser.parse_args( sys.argv[1:] )



def my_rstrip(string, strip):
    if (string.endswith(strip)):
        return string[:-len(strip)]
    return string





alpha_1 = list("ARNDCQEGHILKMFPSTWYV-")
states = len(alpha_1)
alpha_3 = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
           'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','GAP']

aa_1_N = {a:n for n,a in enumerate(alpha_1)}
aa_3_N = {a:n for n,a in enumerate(alpha_3)}
aa_N_1 = {n:a for n,a in enumerate(alpha_1)}
aa_1_3 = {a:b for a,b in zip(alpha_1,alpha_3)}
aa_3_1 = {b:a for a,b in zip(alpha_1,alpha_3)}

def AA_to_N(x):
    # ["ARND"] -> [[0,1,2,3]]
    x = np.array(x);
    if x.ndim == 0: x = x[None]
    return [[aa_1_N.get(a, states-1) for a in y] for y in x]

def N_to_AA(x):
    # [[0,1,2,3]] -> ["ARND"]
    x = np.array(x);
    if x.ndim == 1: x = x[None]
    return ["".join([aa_N_1.get(a,"-") for a in y]) for y in x]


def init_seq_optimize_model():

    hidden_feat_dim = 256
    num_layers = 3

    model = mpnn_util.Struct2Seq(num_letters=21, node_features=hidden_feat_dim, edge_features=hidden_feat_dim, hidden_dim=hidden_feat_dim, num_encoder_layers=num_layers,num_decoder_layers=num_layers, use_mpnn=True, protein_features=args.protein_features, augment_eps=args.augment_eps, k_neighbors=args.num_connections)
    model.to(mpnn_util.device)
    print('Number of model parameters: {}'.format(sum([p.numel() for p in model.parameters()])))
    checkpoint = torch.load(args.checkpoint_path, map_location=torch.device(mpnn_util.device))
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()

    return model



def parse_PDB_biounits(x, atoms=['N','CA','C','O'], chain=None):
    '''
    input:  x = PDB filename
            atoms = atoms to extract (optional)
    output: (length, atoms, coords=(x,y,z)), sequence
    '''
    xyz,seq,min_resn,max_resn = {},{},1e6,-1e6
    for line in x.split("\n".encode("ascii")):
        line = line.decode("utf-8","ignore").rstrip()
  
        if line[:6] == "HETATM" and line[17:17+3] == "MSE":
            line = line.replace("HETATM","ATOM  ")
            line = line.replace("MSE","MET")
    
        if line[:4] == "ATOM":
            ch = line[21:22]
            if ch == chain or chain is None:
                atom = line[12:12+4].strip()
                resi = line[17:17+3]
                resn = line[22:22+5].strip()
                x,y,z = [float(line[i:(i+8)]) for i in [30,38,46]]
        
                if resn[-1].isalpha():
                    resa,resn = resn[-1],int(resn[:-1])-1
                else:
                    resa,resn = "",int(resn)-1
#                 resn = int(resn)
                if resn < min_resn:
                    min_resn = resn
                if resn > max_resn:
                    max_resn = resn
                if resn not in xyz:
                    xyz[resn] = {}
                if resa not in xyz[resn]:
                    xyz[resn][resa] = {}
                if resn not in seq:
                    seq[resn] = {}
                if resa not in seq[resn]:
                    seq[resn][resa] = resi
        
                if atom not in xyz[resn][resa]:
                  xyz[resn][resa][atom] = np.array([x,y,z])
    
    # convert to numpy arrays, fill in missing values
    seq_,xyz_ = [],[]
    try:
        for resn in range(min_resn,max_resn+1):
            if resn in seq:
                for k in sorted(seq[resn]): seq_.append(aa_3_N.get(seq[resn][k],20))
            else: seq_.append(20)
            if resn in xyz:
                for k in sorted(xyz[resn]):
                    for atom in atoms:
                        if atom in xyz[resn][k]: xyz_.append(xyz[resn][k][atom])
                        else: xyz_.append(np.full(3,np.nan))
            else:
                for atom in atoms: xyz_.append(np.full(3,np.nan))
        return np.array(xyz_).reshape(-1,len(atoms),3), N_to_AA(np.array(seq_))
    except TypeError:
        return 'no_chain', 'no_chain'



def generate_seqopt_features( pdb_data, chains ): # multichain
    my_dict = {}
    concat_seq = ''
    concat_N = []
    concat_CA = []
    concat_C = []
    concat_O = []
    concat_mask = []
    coords_dict = {}

    for letter in chains:
        xyz, seq = parse_PDB_biounits(pdb_data, atoms=['N','CA','C','O'], chain=letter)

        concat_seq += seq[0]
        my_dict['seq_chain_'+letter]=seq[0]
        coords_dict_chain = {}
        coords_dict_chain['N_chain_'+letter]=xyz[:,0,:].tolist()
        coords_dict_chain['CA_chain_'+letter]=xyz[:,1,:].tolist()
        coords_dict_chain['C_chain_'+letter]=xyz[:,2,:].tolist()
        coords_dict_chain['O_chain_'+letter]=xyz[:,3,:].tolist()
        my_dict['coords_chain_'+letter]=coords_dict_chain

    my_dict['name']="mpnnv"
    my_dict['num_of_chains'] = len( chains ) 
    my_dict['seq'] = concat_seq

    return my_dict


def sequence_optimize( pdb_data, chains, model ):
    
    t0 = time.time()

    # sequence = get_seq_from_pdb( pdbfile, False )

    feature_dict = generate_seqopt_features( pdb_data, chains )

    #arg_dict = mpnn_util.set_default_args( args.seq_per_struct, decoding_order='random' )
    arg_dict = mpnn_util.set_default_args( 1, decoding_order='random' ) # Setting to 1 since FastRelax protocol doesnt make sense with more cycles
    arg_dict['temperature'] = args.temperature

    masked_chains = [ chains[0] ]
    visible_chains = [ chains[1] ]

    sequences = mpnn_util.generate_sequences( model, feature_dict, arg_dict, masked_chains, visible_chains )
    
    print( f"MPNN generated {len(sequences)} sequences in {int( time.time() - t0 )} seconds" ) 

    return sequences

slash_n = ord("\n")

def process_piped_message(data):
    first_newline = data.find(slash_n)
    if ( first_newline < 0 ):
        return None
    pdb_data = data[first_newline+1:]
    first_line = data[:first_newline].decode("ascii")
    print("MPNNV: ", first_line)

    sp = first_line.split()
    unique_timestamp = sp[0]
    extra_message = ""
    if ( len(sp) > 1 ):
        extra_message = sp[1]

        if ( extra_message == "test" ):
            return (unique_timestamp + " test_complete\n").encode("ascii")

    chains = list(extra_message)

    seq, score = sequence_optimize(pdb_data, chains, mpnn_model )[0]

    return (unique_timestamp + " " + seq + "\n").encode("ascii")



def write_message(f_write, message, add_null=True, max_failures=20):
    written = 0

    max_write = 4096
    if ( add_null ):
        message = message + bytes([0])

    print("MPNNV writing: ", message)
    fail_count = 0
    while ( written < len(message) ):

        we_wrote, written = write_message_one_write(f_write, written, message)

        if ( we_wrote == 0 ):
            time.sleep(0.1)
            fail_count += 1
        else:
            fail_count = 0

        if ( fail_count > max_failures ):
            print("MPNNV could not write to pipe. Quitting")
            sys.exit(1)

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




mpnn_model = init_seq_optimize_model()

if ( args.test_pdb ):

    with open(args.test_pdb, "rb") as f:
        pdb_data = f.read()

    seq, score = sequence_optimize(pdb_data, ['A', 'B'], mpnn_model )[0]
    print(seq)

    sys.exit()


assert(len(args.pipes) == 4)

f_read, f_write, f_close1, f_close2 = args.pipes
os.close(f_close1)
os.close(f_close2)

data_buf = bytes([])

while True:

    try:
        # print("Pre-read")
        data = os.read(f_read, 1024)
        # print("We read", data)
    except OSError as err:
        if err.errno == socket.errno.EAGAIN or err.errno == socket.errno.EWOULDBLOCK:
            # print(err.errno)
            data = None
        else:
            raise  # something else has happened -- better reraise

    null_detected = False

    if ( not data is None ):
        if ( len(data) > 0 ):
            if ( 0 in data ):
                null_detected = True
            data_buf += data

        if ( len(data) == 0 ):
            print("EOF detected")
            break

    while ( null_detected ):
        null_offset = data_buf.find(0)

        this_data = data_buf[:null_offset]
        data_buf = data_buf[null_offset+1:]

        # with open(str(this_data[0:1].decode("ascii")) + "look_at_me.txt", "wb") as f:
        #     f.write(this_data)

        message = process_piped_message(this_data)

        if ( not message is None ):
            write_message(f_write, message, add_null=True)

        null_detected = data_buf.find(0) != -1


os.close(f_write)
os.close(f_read)


print("MPNNV Exit")



