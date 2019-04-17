import sys, string
import spglib
from numpy import *
from numpy.linalg import *

try:
    POSCAR_File = sys.argv[1]
except:
    sys.exit()

########################################################
#########################################################

def isNumber(s):
  try:
    float(s)
    return True
  except ValueError:
    return False

##################################
##atomicNum dictionary

POSCAR = open(POSCAR_File,'r')

title = POSCAR.readline()

scale = map(float,POSCAR.readline().split())
if (len(scale) == 1):
    scale.append(scale[0])    
    scale.append(scale[0])    

aLattice = map(lambda x:x*scale[0],map(float,POSCAR.readline().split()))
bLattice = map(lambda x:x*scale[1],map(float,POSCAR.readline().split()))
cLattice = map(lambda x:x*scale[2],map(float,POSCAR.readline().split()))


i=0
typeAtom=[]
numAtom=[]
while (i==0):
  check=POSCAR.readline().split()
  for mem in check :
    typeAtom.append(mem)
  if(isNumber(check[0])):
    numAtom = map(int,check)
    i=1

  
totNumAtom = reduce(lambda x,y:x+y, numAtom)
sel=0
i=0
while (i==0):
  tmp = POSCAR.readline()
  tmp0=tmp[0].upper()
  if (tmp0=='D' or tmp0=='C'):
    coordType=tmp0
    i=1
  else :
    sel=1
coordList = []
seldy=[]
for member in range(totNumAtom):
    tmtmp=POSCAR.readline().split()
    coordList.append(map(float,tmtmp[0:3]))
    if sel==1 :
        seldy.append(tmtmp[3:6])
    

L=mat([aLattice,bLattice,cLattice])
D=mat(coordList)
if (coordType == 'D'):
  coordList=coordList
else:
  coordList=squeeze(asarray(D*inv(L)))

D=mat(coordList)

symbol_map = {
    "H":1,
    "He":2,
    "Li":3,
    "Be":4,
    "B":5,
    "C":6,
    "N":7,
    "O":8,
    "F":9,
    "Ne":10,
    "Na":11,
    "Mg":12,
    "Al":13,
    "Si":14,
    "P":15,
    "S":16,
    "Cl":17,
    "Ar":18,
    "K":19,
    "Ca":20,
    "Sc":21,
    "Ti":22,
    "V":23,
    "Cr":24,
    "Mn":25,
    "Fe":26,
    "Co":27,
    "Ni":28,
    "Cu":29,
    "Zn":30,
    "Ga":31,
    "Ge":32,
    "As":33,
    "Se":34,
    "Br":35,
    "Kr":36,
    "Rb":37,
    "Sr":38,
    "Y":39,
    "Zr":40,
    "Nb":41,
    "Mo":42,
    "Tc":43,
    "Ru":44,
    "Rh":45,
    "Pd":46,
    "Ag":47,
    "Cd":48,
    "In":49,
    "Sn":50,
    "Sb":51,
    "Te":52,
    "I":53,
    "Xe":54,
    "Cs":55,
    "Ba":56,
    "La":57,
    "Ce":58,
    "Pr":59,
    "Nd":60,
    "Pm":61,
    "Sm":62,
    "Eu":63,
    "Gd":64,
    "Tb":65,
    "Dy":66,
    "Ho":67,
    "Er":68,
    "Tm":69,
    "Yb":70,
    "Lu":71,
    "Hf":72,
    "Ta":73,
    "W":74,
    "Re":75,
    "Os":76,
    "Ir":77,
    "Pt":78,
    "Au":79,
    "Hg":80,
    "Tl":81,
    "Pb":82,
    "Bi":83,
    "Po":84,
    "At":85,
    "Rn":86,
    "Fr":87,
    "Ra":88,
    "Ac":89,
    "Th":90,
    "Pa":91,
    "U":92,
    "Np":93,
    "Pu":94,
    "Am":95,
    "Cm":96,
    "Bk":97,
    "Cf":98,
    "Es":99,
    "Fm":100,
    "Md":101,
    "No":102,
    "Lr":103,
    "Rf":104,
    "Db":105,
    "Sg":106,
    "Bh":107,
    "Hs":108,
    "Mt":109,
    "Ds":110,
    "Rg":111,
    "Cn":112,
    "Uut":113,
    "Uuq":114,
    "Uup":115,
    "Uuh":116,
    "Uus":117,
    "Uuo":118,
    }

Atoms=[]
for i in range(len(numAtom)):
  for j in range(numAtom[i]):
    Atoms.append(symbol_map[typeAtom[i]])

cell = (L,D,Atoms)

#print cell

unit_cell = spglib.refine_cell(cell)
#print unit_cell

prim_cell2 = spglib.find_primitive(cell)
#print prim_cell2
print prim_cell2[0]
print prim_cell2[1]
print prim_cell2[2]

