import argparse
import numpy as np
import re
import warnings
import matplotlib
import matplotlib.pyplot as plt
#import matplotlib.font_manager as font_manager
import json
from matplotlib.collections import LineCollection
from adjustText import adjust_text
#from matplotlib.colors import to_rgba # https://matplotlib.org/2.0.2/api/colors_api.html#matplotlib.colors.to_rgba
# derived from CIAlign
parser = argparse.ArgumentParser(description='Pass a file name.')

parser.add_argument('--f', '--file',  required = True, help = 'input alignment file' )
parser.add_argument('--o', '--output_file', required = True, help = 'output image')
parser.add_argument('--s', '--segments', required = False, help = 'JSON of a hash of notes and MSA coordinate positions of vertical lines')
parser.add_argument('--t', '--title', required = True, help = 'image title')
parser.add_argument('--x', '--xlabel', required = False, help = 'x-axis label', default = 'Amino Acid Residue')
parser.add_argument('--y', '--ylabel', required = False, help = 'y-axis label', default = 'Protein & Species')
args = parser.parse_args()

def base():
    '''
    Returns the hexadecimal values for black and white

    Parameters
    ----------
    None

    Returns
    -------
    dict
        A dictionary containing the hexadecimal values for black and white
    '''

    return {'black': '#000000',
            'white': '#FFFFFF'}

def CBSafe():
    '''
    Returns the hexadecimal values for a colour blind safe colour palette

    Parameters
    ----------
    None

    Returns
    -------
    dict
        A dictionary containing the hexadecimal values for the colours used
        in the CIAlign mini alignments
    '''

    b = base()
    b.update({'yellow_nt': "#c9c433",
              'green_nt': "#56ae6c",
              'red_nt': "#a22c49",
              'blue_nt': "#0038a2",
              'grey_nt': "#6979d3",
              'red_aa': "#a22c49",
              'yellow_aa': "#c9c433",
              'blue_aa':  "#0038a2",
              'orange_aa': "#e57700",
              'midblue_aa': "#589aab",
              'cyan_aa': "#50d3cb",
              'lightgrey_aa': '#eae2ea',
              'green_aa': "#56ae6c",
              'darkgrey_aa': "#888988",
              'purple_aa': '#89236a',
              'paleblue_aa': '#e669ca',
              'peach_aa': "#ffc4a9",
              'tan_aa': "#936e23",
              'remove_insertions': "#9db341",
              'remove_divergent': "#7066bc",
              'crop_ends': '#020545',
              'remove_gap_only': '#f9c1d2',
              'remove_short': "#c85133",
              'crop_divergent': '#ff00d1',
              'user': '#fff4a3'})
    return (b)


def Bright():
    '''
    Returns the hexadecimal values for a colour blind safe colour palette

    Parameters
    ----------
    None

    Returns
    -------
    dict
        A dictionary containing the hexadecimal values for the colours used
        in the CIAlign mini alignments
    '''
    b = base()
    b.update({'yellow_nt': "#ffd500",
              'red_nt': "#64bc3c",
              'green_nt': "#f20707",
              'blue_nt': "#0907f2",
              'grey_nt': "#c7d1d0",
              'red_aa': "#f20707",
              'yellow_aa': "#ffd500",
              'blue_aa':  "#0907f2",
              'orange_aa': "#f4aa03",
              'midblue_aa': "#03b5f4",
              'cyan_aa': "#03f4dd",
              'lightgrey_aa': '#f0f4f4',
              'green_aa': "#64bc3c",
              'darkgrey_aa': "#808080",
              'purple_aa': '#6f0cac',
              'paleblue_aa': '#cde3f8',
              'peach_aa': "#f8e7cd",
              'tan_aa': "#d2a867",
              'remove_insertions': "#9db341",
              'remove_divergent': "#7066bc",
              'crop_ends': '#020545',
              'remove_gaponly': '#f9c1d2',
              'remove_short': "#c85133",
              'crop_divergent': '#ff00d1',
              'user': '#fff4a3'})
    return (b)
try:
    import CIAlign.palettes as palettes
except ImportError:
    import palettes
matplotlib.use('Agg')

def replaceUbyT(arr, rev):
    '''
    Replaces all Us by Ts in the alignment.

    Parameters
    ----------
    arr: np.array
        2D numpu array with the alignment.

    Returns
    -------
    arr: np.array
        2D numpy array of sequences with Ts instead of Us.
    '''
    if not rev:
        arr = np.where(arr == "T", "U", arr)
    else:
        arr = np.where(arr == "U", "T", arr)
    return (arr)

def unAlign(arr):
    '''
    Removes all gaps from the alignment.

    Parameters
    ----------
    arr: np.array
        2D numpu array with the alignment.

    Returns
    -------
    arr: np.array
        2D numpy array of sequences without any gaps.
    '''

    arr = np.where(arr == "-", "", arr)
    return (arr)

def FastaToArray(infile, log=None, outfile_stem=None):
    '''
    Convert an alignment into a numpy array.

    Parameters
    ----------
    infile: string
        path to input alignment file in FASTA format
    log: logging.Logger
        An open log file object

    Returns
    -------
    arr: np.array
        2D numpy array in the same order as fasta_dict where each row
        represents a single column in the alignment and each column a
        single sequence.
    nams: list
        List of sequence names in the same order as in the input file
    '''

    formatErrorMessage = "The MSA file needs to be in FASTA format."
    nams = []
    seqs = []
    nam = ""
    seq = ""
    psl = 0
    nseq = 0
    with open(infile) as input:
        for line in input:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == ">":
                sl = len([s.upper() for s in seq])
                if sl != psl and nseq > 1:
                    print (nseq, sl, psl)
                    raise ValueError ("""
ERROR: The sequences you provided may not be aligned - all the sequences \
are not the same length""")

                psl = sl
                nseq += 1
                seqs.append([s.upper() for s in seq])
                nams.append(nam)
                seq = []
                nam = line.replace(">", "")
            else:
                if len(nams) == 0:
                    if log:
                        log.error(formatErrorMessage)
                    print(formatErrorMessage)
                    exit(1)
                seq += list(line)
    sl = len([s.upper() for s in seq])
    if sl != psl and nseq > 1:
        print (nseq, sl, psl)
        raise ValueError ("""
ERROR: The sequences you provided may not be aligned - all the sequences \
are not the same length""")    
    seqs.append(np.array([s.upper() for s in seq]))
    nams.append(nam)
    arr = np.array(seqs[1:])
    return (arr, nams[1:])


def getPalette(palette='CBS'):
    '''
    Generates a dictionary which assigns a name to each colour using a colour
    blindness safe palette, generated using
    https://medialab.github.io/iwanthue/
    Parameters
    ----------
    palette: str
        The ID of the palette to be used, currently only colour blind safe
        (CBS) is implemented.

    Returns
    -------
    dict
        Dictionary where keys are names of colours and
        values are hexadecimal codes for colours
    '''
    if palette.lower() == 'cbs':
        p = palettes.CBSafe()
    if palette.lower() == 'bright':
        p = palettes.Bright()
    return (p)

def getAAColours(palette='CBS'):
    '''
    Generates a dictionary which assigns a colour to each amino acid.
    Based on the "RasmMol" amino colour scheme and the table here:
        http://acces.ens-lyon.fr/biotic/rastop/help/colour.htm
    approximated using a CB safe palette generated using
    https://medialab.github.io/iwanthue/

    Parameters
    ----------
    pal: str
        A string designating which palette to use, currently only colour blind
        safe (CBS) is implemented.

    Returns
    -------
    dict
        Dictionary where keys are single letter amino acid codes and
        values are hexadecimal codes for colours
    '''
    pal = getPalette(palette=palette)
    return {'D': pal['red_aa'],
            'E': pal['red_aa'],
            'C': pal['yellow_aa'],
            'M': pal['yellow_aa'],
            'K': pal['blue_aa'],
            'R': pal['blue_aa'],
            'S': pal['orange_aa'],
            'T': pal['orange_aa'],
            'F': pal['midblue_aa'],
            'Y': pal['midblue_aa'],
            'N': pal['cyan_aa'],
            'Q': pal['cyan_aa'],
            'G': pal['lightgrey_aa'],
            'L': pal['green_aa'],
            'V': pal['green_aa'],
            'I': pal['green_aa'],
            'A': pal['darkgrey_aa'],
            'W': pal['purple_aa'],
            'H': pal['paleblue_aa'],
            'P': pal['peach_aa'],
            'X': pal['black'],
            '-': pal['white'],
            'B': pal['tan_aa'],
            'Z': pal['tan_aa'],
            'J': pal['tan_aa'],
            '*': pal['white'],
            'U': pal['tan_aa'],
            'O': pal['tan_aa']
            }


def getNtColours(palette='CBS'):
    '''
    Generates a dictionary which assigns a colour to each nucleotide (plus grey
    for "N" and white for "-")
    Parameters
    ----------
    pal: str
        A string designating which palette to use, currently only colour blind
        safe (CBS) is implemented.

    Returns
    -------
    dict
        Dictionary where keys are single letter nucleotide codes and
        values are hexadecimal codes for colours
    '''
    pal = getPalette(palette=palette)
    return {'A': pal['green_nt'],
            'G': pal['yellow_nt'],
            'T': pal['red_nt'],
            'C': pal['blue_nt'],
            'N': pal['grey_nt'],
            "-": pal['white'],
            "U": pal['red_nt'],
            "R": pal['grey_nt'],
            "Y": pal['grey_nt'],
            "S": pal['grey_nt'],
            "W": pal['grey_nt'],
            "K": pal['grey_nt'],
            "M": pal['grey_nt'],
            "B": pal['grey_nt'],
            "D": pal['grey_nt'],
            "H": pal['grey_nt'],
            "V": pal['grey_nt'],
            "X": pal['grey_nt']}

def seqType(arr):
    '''
    Detects if an alignment is of nucleotides or amino acids using pre-built
    dictionarys of amino acid and nucleotide codes.
    Checks if arr contains characters that are not in the dictionary (not
    IUPAC)

    Parameters
    ----------
    arr: np.array
        Numpy array containing the alignment

    Returns
    -------
    str
    'aa' for amino acid and 'nt for nucleotide
    '''
    nt_count = 0
    aa_count = 0
    for seq in arr:
        nucs = set(list(getNtColours().keys()))
        aas = set(list(getAAColours().keys()))
        n = 0
        a = 0
        x = 0
        for s in seq:
            s = s.upper()
            if s in nucs:
                n += 1
            if s in aas:
                a += 1
            if s not in aas and s not in nucs:
                x += 1
        ch = 0
        if n == len(seq):
            nt_count += 1
            ch += 1
        if a == len(seq):
            aa_count += 1
            ch += 1
        if ch == 0:
            print("Unknown nucleotides or amino acids detected.\
                  Please fix your MSA.")
            exit(1)
    if nt_count == len(arr):
        return "nt"
    if aa_count == len(arr):
        return "aa"

def updateNams(nams, removed_seqs):
    '''
    Takes nams, a list of sequence names in the input file, and removed_seqs,
    a set of sequences which have been removed from the file and subtracts
    the removed sequences from nams while maintaining the order.

    Parameters
    ----------
    nams: list
        A list of sequence names in the input file
    removed_seqs: set
        Set of sequence names to remove

    Returns
    -------
    nams2: list
        A list of sequence names which have not been removed
    '''
    nams2 = []
    for nam in nams:
        if nam not in removed_seqs:
            nams2.append(nam)
    return (nams2)


def checkArrLength(arr, log):
    '''
    Checks the shape of the array containing the alignment to ensure that it
    all the sequences haven't been removed by parsing and that all the
    sequences have the same number of columns

    Parameters
    -----------
    arr: np.array
        Numpy array containing the multiple sequence alignment
    log: logging.Logger
        An open log file object
    Returns
    -------
    None
    '''
    emptyAlignmentMessage = """Error: Parsing your alignment with these \
settings has removed all of the sequences."""
    differentLengthsMessage = """Error: The sequences in your alignment are \
not all the same length."""
    if 0 in np.shape(arr):
        log.error(emptyAlignmentMessage)
        print(emptyAlignmentMessage)
        exit(1)
    if len(np.shape(arr)) == 1:
        log.error(differentLengthsMessage)
        print(differentLengthsMessage)
        exit(1)

def arrNumeric(arr, typ, palette='CBS'):
    '''
    Converts the sequence array into a numerical matrix and a colour map
    which matplotlib can interpret as an image (similar to
                                                https://bit.ly/2CIKOEr)
    The rows in the array are inverted so that the output image has the rows
    in the same order as the input alignment.

    Parameters
    ----------
    arr: np.array
        The alignment stored as a numpy array

    typ: str
        Either 'aa' - amino acid - or 'nt' - nucleotide

    palette: str
        Colour palette, CBS or Bright

    Returns
    -------
    arr2: np.array
        The flipped alignment as an array of integers
    cmap: matplotlib.colors.ListedColormap
        A colour map with the colours corresponding to each base
        or amino acid
    '''
    # turn the array upside down
    arr = np.flip(arr, axis=0)
    if typ == 'nt':
        D = getNtColours(palette)
    else:
        D = getAAColours(palette)

    # retrieve the colours for the colour map
    keys = list(D.keys())
    ali_height, ali_width = np.shape(arr)

    # make a dictionary where each integer corresponds to a base or nt
    i = 0
    nD = dict()
    colours = []
    for key in keys:
        if key in arr:
            nD[key] = i
            colours.append(D[key])
            i += 1

    arr2 = np.empty([ali_height, ali_width])

    for x in range(ali_width):
        for y in range(ali_height):
            # numeric version of the alignment array
            arr2[y, x] = nD[arr[y, x]]

    cmap = matplotlib.colors.ListedColormap(colours)
    return (arr2, cmap)
#font = font_manager.FontProperties(style='italic')
def draw_alignment(
		aln_file, # MSA alignment file, probably by clustalo
		height = 3,
		outfile = None,
		segments = None, # will be a JSON string if defined
		title = 'Default Title',
		width = 12,
		xlabel = 'xlabel',
		ylabel = 'Protein & Species'
	):
	if outfile == None:
		raise ValueError('output svg must be named, there is no default')
	arr, names = FastaToArray(aln_file)
	ali_height, ali_width = np.shape(arr)
	arr2, cm = arrNumeric(arr, seqType(arr), 'CBS')
	lineweight_h = 10 / ali_height
	lineweight_v = 10 / ali_width

	f = plt.figure(figsize=(width, height))
	a = f.add_subplot(1, 1, 1)
	
	a.set_xlim(-0.5, ali_width)
	a.spines['right'].set_visible(False)
	a.spines['top'].set_visible(False)
	a.spines['left'].set_visible(False)
	f.suptitle(title, fontsize=11)#, y=0.92)
	# display it on the axis
	a.imshow(arr2, cmap=cm, aspect='auto', interpolation='nearest')
	list_yticks = list(range(0,len(names)))
	list_yticks.reverse()
	a.set_yticks(list_yticks)
	a.set_yticklabels(names, fontsize = 10)
	a.set_xlabel(xlabel)
	a.set_ylabel(ylabel)
	# add segments if they're defined
	if segments != None: # add vertical lines (line segments) to describe active sites
		active_sites = json.loads(segments)
		lines = []
		colors = []
		texts = [a.text(active_sites[name], -1, name) for name in active_sites]
		for name in active_sites:
			lines.append([(active_sites[name], len(arr)), (active_sites[name], 0)])
			colors.append([0.0, 0.0, 0.0, 1.0]) # rgb for black color
		lc = LineCollection(lines, colors = colors, linewidth = 3, alpha = 0.5, linestyle = 'dashed')
		a.add_collection(lc)
		adjust_text(texts, only_move={"text": "x", "static": "x", "explode": "x", "pull": "x"})
		#a.autoscale()
	f.savefig(outfile, bbox_inches='tight', pad_inches = 0.1, metadata={'Creator': 'made/written by ' + __file__ })
	print('wrote ' + outfile)
	
draw_alignment(
	aln_file = args.f,
	title    = args.t,
	outfile  = args.o,
	segments = args.s,
	xlabel   = args.x,
	ylabel   = args.y
)
