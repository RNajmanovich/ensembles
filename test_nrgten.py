import os
from pathlib import Path
from inspect import currentframe, getframeinfo
import warnings
import argparse
from nrgten.encom import ENCoM # needed for block 1
from Bio import SeqIO # needed for block 2
# from Bio import Align # needed for block 3
from Bio.PDB import * # block 3
from modeller import * # needed for block 4
from modeller.automodel import * # block 5

# this is a test
def main():
    parse_initializer = argparse.ArgumentParser()
    parse_initializer.add_argument("-f", "--filein", help = "Input PDB structure filename", type=str, required=True)
    # parse_initializer.add_argument("-d","--directory", help = "Directory for output", type=str, const = '.', default = '.') not implemented yet
    parse_initializer.add_argument("-w", "--warnings", help = "Do not silence code warnings", action='store_false')
    args = parse_initializer.parse_args()
    
    if args.warnings:
        warnings.simplefilter('ignore')


    filein_code = args.filein.split(".")[0]
    # print("The input PDB filename is: %s " % args.filein)
    # print("The input PDB coode is: %s " % filein_code)


    # 1: Runs ENCoM on the input structure to generate a conformational ensemble
    step_val = 1.0
    max_disp_val = 2.0
    ensemble_base = str(step_val) + "_" + str(max_disp_val)
    ensemble_file = ensemble_base + ".pdb"
    model = ENCoM(args.filein)
    model.compute_bfactors()
    model.write_dynamical_signature(ensemble_base + "_signature.txt")
    model.build_conf_ensemble([7,8,9], ensemble_file, step=step_val, max_displacement=max_disp_val)

    # frameinfo = getframeinfo(currentframe())
    # print("paused at: ", frameinfo.filename, frameinfo.lineno)
    # input()

    # 2: gets sequence from a PDB file atom records.
    for record in SeqIO.parse(args.filein, "pdb-atom"):
        print(record.seq)

    # 5: opens the ensemble_file, and for each model in the ensemble, writes a PDB file and uses it 
    #    as template for Modeller

    maxmodels = 1

    p = PDBParser()
    structure = p.get_structure("ensb",ensemble_file)
    total_states = len([_ for _ in structure.get_models()])
    # print(str(total_states))
    # frameinfo = getframeinfo(currentframe())
    # print("paused at: ", frameinfo.filename, frameinfo.lineno)
    # input()

    for encom_state in structure:
        # this block is for testing purposes
        # if encom_state.id > 5:
        #     continue
        

        # creates the encom_*.pdb files each containing a single model of the ensemble created by encom
        io = PDBIO()
        io.set_structure(encom_state)
        print("ID", encom_state.id, "/" + str(total_states) + "\n")
        encom_base = "encom_" + str(encom_state.id)
        encom_fileout_pdb = encom_base + ".pdb"
        io.save(encom_fileout_pdb, preserve_atom_numbering = True)

        # frameinfo = getframeinfo(currentframe())
        # print("paused at: ", frameinfo.filename, frameinfo.lineno)
        # input()

        # 4: Creates the PIR format alignment required for modeller. It takes the sequence from the atom part of the 
        #    original PDB file used to create the ensemble.
        filepir = encom_base + ".pir"
        pir = open(filepir, "w")
        pir.write(">P1;" + encom_base + "\n")
        pir.write("sequence:" + encom_base + ":::::::0.00: 0.00\n") # this line has to be improved to write correct info
        pir.write(str(record.seq) + "*")
        pir.close()

        # frameinfo = getframeinfo(currentframe())
        # print("paused at: ", frameinfo.filename, frameinfo.lineno)
        # input()

        env = Environ()
        aln = Alignment(env)
        mdl = Model(env, file=filein_code)
        aln.append_model(mdl, align_codes=encom_base, atom_files=encom_fileout_pdb)
        aln.append(file=filepir, align_codes=encom_base)
        aln.align2d(max_gap_length=50)
        aln.write(file='tmpl_modl.ali', alignment_format='PIR')
        model_base = "encom_" + str(encom_state.id)
    
        # frameinfo = getframeinfo(currentframe())
        # print("paused at: ", frameinfo.filename, frameinfo.lineno)
        # input()

        env = Environ()
        a = AutoModel(env, alnfile='tmpl_modl.ali', knowns=model_base, sequence=encom_base, assess_methods=(assess.DOPE))
        #a.very_fast()
        #a = AutoModel(env, alnfile='tmpl_modl.ali', knowns='modl', sequence='tmpl')
        a.starting_model = 1
        a.ending_model = maxmodels
        a.make()
        for i in range(maxmodels):
            # print(str(i+1))
            modellername = encom_base + ".B9999000" + str(i+1) + ".pdb"
            newmdlname = "model_" + str(encom_state.id) + "_" + str(i) + ".pdb"
            os.rename(modellername,newmdlname)
        for p in Path(".").glob(encom_base + '.*'):
            if str(p)  != encom_fileout_pdb:
                # print("will delete:" + str(p))
                p.unlink()
        
            # frameinfo = getframeinfo(currentframe())
            # print("paused at: ", frameinfo.filename, frameinfo.lineno)
            # input()
        if os.path.exists('tmpl_modl.ali'):
            os.remove('tmpl_modl.ali')

if __name__ == '__main__':
    main()