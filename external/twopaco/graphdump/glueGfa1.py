import sys
from Bio.Seq import Seq

#An utility that takes a GFA generated by TwoPaCo as input and glues it back into
#original genome(s)

def rev_comp(seq):
	rc = Seq(seq).reverse_complement()
	return str(rc)

if len(sys.argv) != 2:
	print "Usage: glue.py <GFA file>"
	exit()

segment = dict()

def spell_segment(sid):
	return segment[sid] if sid > 0 else rev_comp(segment[-sid])

def spell_path(path, k):
	chr_body = [spell_segment(path[0])]
	spell = [spell_segment(path[0])]
	for i in xrange(1, len(path)):
		sid = path[i]
		s = spell_segment(sid)
		ps = spell[-1]
		if ps[-k:] != s[:k]:
			print path[i], path[i - 1], path
			print spell, s
			sys.exit("Improper segment overlap!")
		spell.append(s)
		chr_body.append(s[k:])
	print ''.join(chr_body)	

handle = open(sys.argv[1])
for line in handle:
	line = line.split()
	if len(line) > 0:
		if line[0] == 'S':
			sid, body = line[1], line[2]
			if sid not in segment:
				if body != '*':
					segment[int(sid)] = body
			else:
				print line
				sys.exit("Segment duplicate!")
		if line[0] == 'P':
			k = int(line[3].split(',')[0][:-1])
			path = [int(x) for x in line[2].split(',')]
			spell_path(path, k)
handle.close()