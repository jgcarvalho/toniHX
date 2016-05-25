#!/usr/bin/env python
import sys,os,math,re
import argparse
from glob import glob

beta_c = 0.35
beta_h = 2.0
hcut = 2.4 # Angstrom
heavycut = 6.5 # Angstrom
excl = 2
hxtarg = -1. # maxwell demon
hxforce = 1000.

def dist(i,j):
    dx=i[2]-j[2]
    dy=i[3]-j[3]
    dz=i[4]-j[4]
    return math.sqrt(dx**2+dy**2+dz**2);

def aproc(x):
    res = int(x[22:26])
    name = x[17:20]
    chain = x[21]
    idx = int(x[6:11])
    X = float(x[30:38])
    Y = float(x[38:46])
    Z = float(x[46:54])
    return (res,idx,X,Y,Z, chain, name)

def test(x):
    print x

def process(fn):
    atoms = filter(lambda x: x[0:6] =="ATOM  ", open(fn).readlines())
    #hxres = map(lambda x: int(x.split()[0]), open(sys.argv[2]).readlines())
    #logp = map(lambda x: float(x.split()[1]), open(sys.argv[2]).readlines())


    N = map(aproc,filter(lambda x: x.split()[2]=="N", atoms));
    H = map(aproc,filter(lambda x: x.split()[2] =="H", atoms));
    O = map(aproc,filter(lambda x: x.split()[2]=="O", atoms));
    # correction to exclude hydrogen atoms
    #heavy = map(test,filter(lambda x: x[13] != "H" and x.split()[2][0] != "H", atoms));
    heavy = map(aproc,filter(lambda x: x[13] != "H" and x.split()[2][0] != "H", atoms));
    hbonds = {}

    hproc = {}
    dist_ho = {}


    for h in H:
        res = (h[5], h[0], h[6])
        # hbonds[res] = [h[1]]
        if res not in hbonds.keys():
            hbonds[res] = 0
        if res not in hproc.keys():
            hproc[res] = 0.
        for o in O:
            if dist(h,o) < hcut:
                # hbonds[res].append(o[1])
                hbonds[res] += 1
                hproc[res] += beta_h/(1.+math.exp(10.*(dist(h,o)-hcut)))
                if res in dist_ho.keys():
                    dist_ho[res].append(dist(h,o))
                else:
                    dist_ho[res] = [dist(h,o)]

    # for res in dist_ho.keys():
    #     print res, dist_ho[res]

    heavycontacts = {}
    dist_nhv = {}

    # heavy atom contacts
    for n in N:
        res = (n[5], n[0], n[6])
        # heavycontacts[res] = [n[1]]
        if res not in heavycontacts.keys():
            heavycontacts[res] = 0
        if res not in hproc.keys():
            hproc[res] = 0.
        for hv in heavy:
            if dist(n,hv) < heavycut and abs(n[0]-hv[0])>=excl:
                # heavycontacts[res].append(hv[1])
                heavycontacts[res] += 1
                hproc[res] += beta_c/(1.+math.exp(5.*(dist(n,hv)-heavycut)))
                if res in dist_nhv.keys():
                    dist_nhv[res].append(dist(n,hv))
                    pass
                else:
                    dist_nhv[res] = [dist(n,hv)]

    # print "heavy"
    # for res in dist_nhv.keys():
    #     print res, dist_nhv[res]
    #print hbonds
    #print heavycontacts

    residues = hbonds.keys()
    residues.sort()
    return (fn,residues, hbonds, heavycontacts, hproc, dist_ho, dist_nhv)

# write the data (d) to a csv file (outfile)
def report(d, outfile):
    f = open(outfile, 'w')
    if ".csv" in outfile:
        fdho = open(outfile.replace(".csv", "_dist_ho.csv"), 'w')
        fnhv = open(outfile.replace(".csv", "_dist_nhv.csv"), 'w')
    else:
        fdho = open(outfile + "_dist_ho", 'w')
        fnhv = open(outfile + "_dist_nhv", 'w')


    for col in range(0,len(d)):
        fn = os.path.basename(d[col][0])
        if col == 0:
            f.write("ResNumber, ResName, Chain, nHB_{}, nC_{}, hproc1_{}, hproc2_{}".format(fn,fn,fn,fn))
            fdho.write("ResNumber, ResName, Chain, dist_HO_{}".format(fn))
            fnhv.write("ResNumber, ResName, Chain, dist_NHeavy_{}".format(fn))
            if len(d) == 1:
                f.write("\n")
                fdho.write("\n")
                fnhv.write("\n")
        elif col == len(d)-1:
            f.write(", nHB_{}, nC_{}, hproc1_{}, hproc2_{}\n".format(fn,fn,fn,fn))
            fdho.write(", dist_HO_{}\n".format(fn))
            fnhv.write(", dist_NHeavy_{}\n".format(fn))
        else:
            f.write(", nHB_{}, nC_{}, hproc1_{}, hproc2_{}".format(fn,fn,fn,fn))
            fdho.write(", dist_HO_{}".format(fn))
            fnhv.write(", dist_NHeavy_{}".format(fn))

    for res in d[0][1]:
        for col in range(0,len(d)):
            # nhb = len(d[col][2][res]) - 1
            nhb = d[col][2][res]
            # nc = len(d[col][3][res]) - 1
            nc = d[col][3][res]
            hp = d[col][4][res]
            if res in d[col][5].keys():
                distHO = d[col][5][res]
            else:
                distHO = []
            if res in d[col][6].keys():
                distNHv = d[col][6][res]
            else:
                distNHv = []

            if col == 0 :
                f.write("{}, {}, {}, {}, {}, {}, {}".format(res[1],res[2],res[0],nhb,nc,nhb*beta_h + nc*beta_c,hp))
                fdho.write("{}, {}, {}, \"{}\"".format(res[1],res[2],res[0],distHO))
                fnhv.write("{}, {}, {}, \"{}\"".format(res[1],res[2],res[0],distNHv))
                if len(d) == 1:
                    f.write("\n")
                    fdho.write("\n")
                    fnhv.write("\n")
            elif col == len(d)-1:
                f.write(", {}, {}, {}, {}\n".format(nhb,nc,nhb*beta_h + nc*beta_c,hp))
                fdho.write(", \"{}\"\n".format(distHO))
                fnhv.write(", \"{}\"\n".format(distNHv))
            else:
                f.write(", {}, {}, {}, {}".format(nhb,nc,nhb*beta_h + nc*beta_c,hp))
                fdho.write(", \"{}\"".format(distHO))
                fnhv.write(", \"{}\"".format(distNHv))
    f.close()
    fdho.close()
    fnhv.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="PDB files directory")
    parser.add_argument("output", help="CSV output file")
    args = parser.parse_args()

    files = glob(args.directory+"/*.pdb")
    outfile = args.output

    all_pdbs = []
    # print files
    for i in files:
        data = process(i)
        all_pdbs.append(data)
    print "Number of PDB files in \"{}\": {}".format(args.directory, len(all_pdbs))
    report(all_pdbs, outfile)
