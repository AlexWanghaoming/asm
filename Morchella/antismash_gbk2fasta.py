from Bio import SeqIO,Seq
import sys,re

def gbk2fasta(fhs):
	for fh in fhs:
		# for record in SeqIO.parse("/Users/alexwang/0data/oyr/sm/contig_4.region003.gbk", format="genbank"):
		for record in SeqIO.parse(fh, format="genbank"):
			for seq in record.features:
				try:
					if seq.qualifiers['gene_kind'] == ["biosynthetic"]:
						aa_id = seq.qualifiers['gene_id'][0]
						aa_seq = seq.qualifiers['translation'][0]
						func = seq.qualifiers['gene_functions'][0]
						capture = re.findall(r"biosynthetic \(rule-based-clusters\) (.*):.*",func)[0]
						# s = SeqIO.SeqRecord(Seq.Seq(aa_seq), id=aa_id, description="")
						# SeqIO.write(s, "asd", format="fasta")
						sys.stdout.write(">{0};{1}\n{2}\n".format(aa_id, capture, aa_seq))
				except KeyError:
					pass

if __name__ == '__main__':
	# gbk2fasta()
	fhs = sys.argv
	gbk2fasta(fhs)
