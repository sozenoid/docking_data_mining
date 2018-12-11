###
### some scripts to extract the data from AutoDock Vina logs
###


def create_dictionary_smi_pbnbr(listoflogs, pubchem_smi_database_file, outf):
    """
    :param listoflogs: takes a list of files named "./83527306.pdbqt-LOG" where the number is a puchem number and the file is formatted as:
    #################################################################
    # If you used AutoDock Vina in your work, please cite:          #
    #                                                               #
    # O. Trott, A. J. Olson,                                        #
    # AutoDock Vina: improving the speed and accuracy of docking    #
    # with a new scoring function, efficient optimization and       #
    # multithreading, Journal of Computational Chemistry 31 (2010)  #
    # 455-461                                                       #
    #                                                               #
    # DOI 10.1002/jcc.21334                                         #
    #                                                               #
    # Please see http://vina.scripps.edu for more information.      #
    #################################################################

    Reading input ... done.
    Setting up the scoring function ... done.
    Analyzing the binding site ... done.
    Using random seed: 42
    Performing search ... done.
    Refining results ... done.

    mode |   affinity | dist from best mode
         | (kcal/mol) | rmsd l.b.| rmsd u.b.
    -----+------------+----------+----------
       1         -2.9      0.000      0.000
       2         -2.9      9.742     13.516
       3         -2.9      9.725     13.596
       4         -2.9      9.699     13.518
       5         -2.9      2.434      2.610
       6         -2.8      9.733     13.555
       7         -2.8      9.568     13.178
       8         -2.8      9.635     13.470
       9         -2.7      2.702      3.553
    Writing output ... done.

    :param pubchem_smi_database_file: takes a single file that matches pubchem numbers to smiles as:
    O(C(C[N+](C)(C)C)CC(=O)[O-])C(=O)C	1
    :param outf: name of the output file for the dictionary
    :return: a dictionary where the keys are pubchem numbers for the files in listoflogs and the values are the associated smiles
    """

    keys = set()
    for f in listoflogs:
        try:
            pbnum = f.split("/")[-1].split(".")[0]
            keys.add(int(pbnum))
        except:
            print f
    dic = dict.fromkeys(keys)

    with open(pubchem_smi_database_file, 'rb') as r:
        for line in r:
            smiles, pbnum = line.strip().split('\t')
            pbnum = int(pbnum)
            if pbnum in dic:
                dic[pbnum] = smiles

    with open(outf, 'wb') as w:
        cPickle.dump(dic, w)


if __name__ == "__main__":
    import glob
    import cPickle
    create_dictionary_smi_pbnbr(glob.glob("/home/macenrola/Documents/docked_for_data_analysis/docked_logs/*LOG"),
                                "/home/macenrola/Documents/docked_for_data_analysis/pubchem_smis",
                                "/home/macenrola/Documents/docked_for_data_analysis/500k_docked_smidic")