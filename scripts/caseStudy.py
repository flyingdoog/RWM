ac = []
pc = []

with open('./casestudy.txt') as fin:
	line = fin.readline()
	ac = [int(s) for s in line.strip().split(' ')]
	line = fin.readline()
	pc = [int(s) for s in line.strip().split(' ')]
print(ac)
print(pc)

acset = set(ac)
pcset = set(pc)

names = []

with open('C:\\Users\\dul262\\Desktop\\MultiNetwork\\DBLP_Citation_2014_May\\original\\a2id.txt') as fin:
	for line in fin:
		name = line.strip().split('\t')
		names.append(name)

papers = []
p2id = {}
with open('C:\\Users\\dul262\\Desktop\\MultiNetwork\\DBLP_Citation_2014_May\\original\\p2id.txt') as fin:
	for line in fin:
		name = line.strip().split('\t')[0]
		p2id[name]=len(papers)
		papers.append(name)


venues = []
titles = []
years = []
authors = []


title = ''
year = ''
venue = ''
pid = 0

with open('C:\\Users\\dul262\\Desktop\\MultiNetwork\\DBLP_Citation_2014_May\\original\\publications.txt', encoding='utf-8') as fin:
	for line in fin:
		if line.startswith('#*'):
			title = line.strip()[2:]
		elif line.startswith('#@'):
			author = line.strip()[2:]
		elif line.startswith('#t'):
			year = line.strip()[2:]
		elif line.startswith('#c'):
			venue = line.strip()[2:]
		elif line.startswith('#index'):
			index = line.strip()[6:]
			pid = p2id[index]
		elif line.startswith('#!'):
			if pid in pcset:
				titles.append(title)
				authors.append(author)
				venues.append(venue)
print(venues)
print(authors)
print(titles)

for aid in ac:
	print(names[aid])

