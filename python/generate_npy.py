from mdance.inputs.preprocess import gen_traj_numpy
import sys
import argparse
import numpy as np
from pathlib import Path

def parser(args : list[str]):
    p = argparse.ArgumentParser("Generate Numpy")
    p.add_argument('-top', '--topology', type=str, metavar='[PATH]', 
                   required=True, help="Input topology file", dest='top')
    p.add_argument('-traj', '--trajectory', type=str, metavar='[PATH]',
                   required=True, help="Input trajectory file", dest='traj')
    p.add_argument('-out', '--output', type=str, default='output', metavar='[OUTPUT NAME]',
                   help='Output file path without file extension', dest='output')
    p.add_argument('-sel', '--selection', nargs='*', default='', dest='sel',
                   help='Atom selection string (e.g. resid 3 to 12 and name N CA C O H)')
    return p.parse_known_args(args)

def main():
    args = parser(sys.argv[1:])[0]
    input_top = args.top
    input_traj = args.traj
    output_base_name = args.output
    atomSelection = " ".join(args.sel).strip('\'\"') if type(args.sel) == list else args.sel

    traj_numpy = gen_traj_numpy(input_top, input_traj, atomSelection)

    Path(output_base_name).parent.mkdir(parents=True, exist_ok=True)
    output_name = output_base_name + '.npy'
    np.save(output_name, traj_numpy)

if __name__ == "__main__":
    main()