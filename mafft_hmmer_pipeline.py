
# Iterates through multi-fasta files, aligns them with MAFFT, creates hmm profiles and compares them to a genome

import subprocess as sp
from optparse import OptionParser
import os
import glob


parser = OptionParser()

parser.add_option('-a', '--mafft', dest='alignment',
                  help='Align multifasta files?',
                  metavar='<Yes/ No>')
parser.add_option('-m', '--hpr', dest='hmmpr',
                  help='Do you want to do a HMM profiles of the alignment?.',
                  metavar='<Yes/ No>')
parser.add_option('-n', '--cpu', dest='cpu',
                  help='Number of CPUs to use',
                  metavar='<n>')
parser.add_option('-k', '--hmms', dest='hmmsearch',
                  help='Search profiles coincidences in databases?',
                  metavar='<Yes/ No>')
parser.add_option('-E', '--Eval', dest='Evalue',
                  help=' E-value for the HMM search. (mandatory if -m Yes',
                  metavar='<ne-m>')
(options, args) = parser.parse_args()


# Function definition


# Iterates through multifasta files and align them with MAFFT. It requires multifasta files to be in
# a folder named 'multifasta'.

def align_mafft(multifasta_file):

    if not os.path.exists('protein_alignments'):

        os.makedirs('protein_alignments')
        pass

    if not os.path.exists('dna_alignments'):

        os.makedirs('dna_alignments')
        pass

    filename = multifasta_file.split('/', 2)[1].split('.fasta', 2)[0]

    pmafft = 'mafft --thread -1 --auto --globalpair multifasta/%s.fasta > protein_alignments/%s.fasta' % (filename, filename)
    dmafft = 'mafft --thread -1 --auto --globalpair multifasta/%s.fasta > dna_alignments/%s.fasta' % (filename, filename)

    if multifasta_file.startswith('multifasta/p'):

        sp.call(pmafft, shell=True)

    else:

        sp.call(dmafft, shell=True)

# Build a HMM profile for each alignment

def make_hmm(alnfile):

    if not os.path.exists('profiles'):

        os.makedirs('profiles')
        pass

    filename = alnfile.split('/', 2)[1].split('.fasta', 2)[0]

    if options.cpu is not None:

        hmmcmd = ('hmmbuild --cpu %s profiles/%s.hmm ' + '%s') % (options.cpu, filename, alnfile)

    else:

        hmmcmd = ('hmmbuild profiles/%s.hmm ' + '%s') % (filename, alnfile)
        pass

    sp.call(hmmcmd, shell=True)
    pass

# Search coincidences in a database using previously created HMM profiles as a query.
# It requires assembled transcriptomes names to start with 'SRR' to search nucleotide coincidences.
# It requires genome protein files names to start with 'GCA' to search protein coincidences.

def search_hmm(profile, seqdb):

    if not os.path.exists('searches'):

        os.makedirs('searches')
        pass

    seqdb_name = seqdb.split('/')[1].split('.f')[0]
    print('Database: ' + seqdb_name)

    if not os.path.exists('searches/%s' % seqdb_name):

        os.makedirs('searches/%s' % seqdb_name)
        pass

    filename = seqdb_name + '_' + profile.split('/', 2)[1].split('.hmm', 2)[0]
    print(filename)
    hpr_cmd = 'hmmpress %s' % profile

    sp.call(hpr_cmd, shell=True)

    if options.cpu is not None:

        hmmcmd = 'hmmsearch --cpu %s -E %s -o searches/%s/%s_hmmsearch.txt %s %s ' % \
                 (options.cpu, options.Evalue, seqdb_name, filename, profile, seqdb)

    else:

        hmmcmd = 'hmmsearch -E %s -o s%s_hmmsearch.txt %s %s' % \
                 (options.Evalue, seqdb_name, filename, profile, seqdb)

        pass

    if filename.startswith('SRR') & filename.__contains__('1_d'):

        sp.call(hmmcmd, shell=True)
        sp.call('rm profiles/*.hmm.*', shell=True)

    elif filename.startswith('GCA') & filename.__contains__('n_p'):

        sp.call(hmmcmd, shell=True)
        sp.call('rm profiles/*.hmm.*', shell=True)


    pass

# Function calling

if options.alignment == 'Yes' or options.alignment == 'yes' or options.alignment == 'Y':

    multifastas = glob.glob('multifasta/*')

    for multifasta in multifastas:

        align_mafft(multifasta)
        pass

    pass

if options.hmmpr == 'Yes' or options.hmmpr == 'yes' or options.hmmpr == 'Y':

    files = glob.glob('protein_alignments/*')
    files.extend(glob.glob('dna_alignments/*'))

    for file in files:

        make_hmm(file)
        pass

    pass

if options.hmmsearch == 'Yes' or options.hmmsearch == 'yes' or options.hmmsearch == 'Y':

    profiles = glob.glob('profiles/*.hmm')
    DBs = glob.glob('DBs/*.fa*')

    for profile in profiles:

        for DB in DBs:

            search_hmm(profile, DB)
            pass

        pass

    pass

else:
    print('Error')
    pass
