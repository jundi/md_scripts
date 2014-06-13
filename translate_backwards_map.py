import sys

# input/output files
#--------------------------
oldmapfile = 'popc.charmm36.map'
newmapfile = 'popc.opls.map'
itpfile = 'popc-opls.itp'
transfile = 'popc-charm2opls.trans'


# load translation
#--------------------------
transstream = open(transfile, 'r')
newtoold = {}
oldtonew = {}

for line in transstream:

  if line.startswith(';'):
    continue

  columns = line.split()
  if len(columns) < 2:
    continue

  newtoold[columns[1]] = columns[0]
  oldtonew[columns[0]] = columns[1]


transstream.close()



# load atom order from itp
#--------------------------
itpstream = open(itpfile, 'r')
order = []

# skip stuf  before [ atoms ] section
while itpstream.readline().find('[ atoms ]') < 0:
  continue

# skip header
itpstream.readline()

# write [ atoms ] section to tempitp
while True:
  line = itpstream.readline()
  if line.startswith('['):
    break
  if line.startswith(';'):
    continue

  columns = line.split()
  if len(columns) < 5:
    continue

  atom = columns[4]
  order.append(atom)

itpstream.close()



# load cg-beads from old map file
#--------------------------
oldmapstream = open(oldmapfile, 'r')
beads = {}

# skip stuf  before [ atoms ] section
while oldmapstream.readline().find('[ atoms ]') < 0:
  continue

# save beads
while True:
  line = oldmapstream.readline()
  if line.startswith('['):
    break
  if line.startswith(';'):
    continue

  columns = line.split()
  if len(columns) < 3:
    continue

  atom = columns[1]
  beadset = columns[2:]
  beads[atom] = beadset

oldmapstream.close()




# write new file
#--------------------------
oldmapstream = open(oldmapfile, 'r')
newmapstream = open(newmapfile, 'w')

# Until [ atoms ] section, just copy everything
while True:
  line = oldmapstream.readline()
  newmapstream.write(line)
  if line.find('atoms') > 0:
    break
  if len(line) < 1:
    print '[ atoms ] section not found.'
    sys.exit(0)


# write [ atoms ] section
num = 1
for atom in order:
   
  oldatom = newtoold[atom]
  beadstr = ' '.join(beads[oldatom])
  atomstr = "{0:5}{1:6}".format(str(num), atom)
  newline = atomstr+ '   ' + beadstr + '\n'
  newmapstream.write(newline)
  num = num + 1


# find next section after [ atoms ]
while True:
  line = oldmapstream.readline()
  if line.startswith(';') or line.startswith('[') or line.isspace():
    newmapstream.write(line)
    break

# Switch atom names in the last sections
while True:

  line = oldmapstream.readline()

  if line.startswith(';') or line.startswith('[') or line.isspace():
    newmapstream.write(line)
    continue

  if len(line) < 1:
    print 'end of file'
    break

  newline=''
  oldatoms = line.split()
  for a in oldatoms:
    newline = newline + oldtonew[a] + ' '

  newmapstream.write(newline + '\n')


oldmapstream.close()
newmapstream.close()
