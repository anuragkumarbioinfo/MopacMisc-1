import sys
inaux = sys.argv[1] 
inref = sys.argv[2]

print "PHASE I: AUX CODE INITIATED"

RADIUS = {'Ni': '1.63', 'Tl': '1.96', 'As': '1.85', 'O': '1.52', 'Au': '1.66', 
          'Ar': '1.88', 'He': '1.4', 'Pd': '1.63', 'Se': '1.9', 'Zn': '1.39', 
          'Si': '2.1', 'Ga': '1.87', 'Sn': '2.17', 'Xe': '2.16', 'In': '1.93', 
          'Ne': '1.54', 'Na': '2.27', 'N': '1.55', 'Te': '2.06', 'Pt': '1.75', 
          'U': '1.86', 'Cl': '1.75', 'Li': '1.82', 'Br': '1.85', 'H': '1.2', 
          'F': '1.47', 'Pb': '2.02', 'S': '1.8', 'K': '2.75', 'Cu': '1.4', 
          'Kr': '2.02', 'I': '1.98', 'P': '1.8', 'Ag': '1.72', 'Mg': '1.73', 
          'C': '1.7', 'Cd': '1.58', 'Hg': '1.55'}

def RADIUS2(x):
    if x not in RADIUS: 
        print("No Van der waal radius found for Element: "+x)
        return (2.0)
    else: 
        return RADIUS[x]

ELE = "ATOM_EL"
CHR = "ATOM_CHARGES"
OPT = "ATOM_X_"
OPTE = "HEAT_OF_FORM"
e,c,o = False, False, False

print "PHASE II: Reading Reference PDB FILE"


#PQR = "ATOM    {:3}  {:4}{:3.7}  {:3}   {:8.3f}{:8.3f}{:8.3f}{:8.2f}{:7.2f}"
#TODO: Fix this. Some problem in Atom number spacing.
PQR =  "ATOM    {:3} {:4} {:3.7}   {:3}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}"
#ATOM     28  C   LEU B 466     118.113 112.495  46.722  1.00  0.00           C  
refPDB = [(i[17:20],i[23:26]) for i in open(inref).readlines() if i.startswith("ATOM") or i.startswith("HETATM")]

HOLD_VALUES = {"XYZ" : {}, 
               "CHARGE" : [],
               "ELEMENT" : []}

if "old" in sys.argv:
    OPT = "ATOM_X:ANGSTROMS"
    print("Warning: Switching to old coordinates (non-optimized)")

opt = 0

print "PHASE III: Reading AUX FILE"

with open(inaux) as a:
    for num, lines in enumerate(a):
        #if "NOOPT" in lines:
        #    OPT = "ATOM_X:ANGSTROMS"
        if ELE in lines:
            e = True
        elif CHR in lines:
            c = True
        elif OPTE in lines:
            o = False
        elif OPT in lines:
            o = True
            opt += 1
            if opt not in HOLD_VALUES["XYZ"]: HOLD_VALUES["XYZ"][opt] = []
        elif "=" in lines:
            e,c,o = False, False, False
        elif e:
            HOLD_VALUES["ELEMENT"]+= lines.strip().split()
        elif c:
            HOLD_VALUES["CHARGE"]+=list(map(float, lines.strip().split()))
        elif o:
            try:
                HOLD_VALUES["XYZ"][opt].append(list(map(float,lines.strip().split())))
            except:
                print lines
                pass
        #print num

#print HOLD_VALUES["XYZ"].keys()

#raise SystemExit
print "PHASE IV: Writing PQR FILE"

with open(inaux.replace(".aux",'.out.pqr'), "w") as output:
        output.write("MODEL    1\n")
        for model in range(1, len(HOLD_VALUES["XYZ"].keys()) + 1):
            for atom in range(len(HOLD_VALUES["XYZ"][1])):
                element = HOLD_VALUES["ELEMENT"][atom]
                try:
                    rn, rid = refPDB[atom]
                    output.write(PQR.format(atom +1,
                                             element+str(atom +1),
                                             rn, #RESIDUE NAME
                                             int(rid), #RESIDUE ID
                                             HOLD_VALUES["XYZ"][model][atom][0],
                                            HOLD_VALUES["XYZ"][model][atom][1],
                                             HOLD_VALUES["XYZ"][model][atom][2],
                                            float(HOLD_VALUES["CHARGE"][atom]),
                                             float(RADIUS2(element)) )+ "\n")
                except IndexError:
                    output.write(PQR.format(atom +1, element+str(atom +1), "AUX", 1, HOLD_VALUES["XYZ"][model][atom][0],HOLD_VALUES["XYZ"][model][atom][1], HOLD_VALUES["XYZ"][model][atom][2], 
                                            float(0.0), float(RADIUS2(element)) )+ "\n")
                
                
            output.write("ENDMDL\nMODEL    {}\n".format(model+1))

print "PHASE V: Finished!"
