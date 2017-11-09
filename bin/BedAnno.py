#!/usr/bin/env python
# coding: utf-8
from operator import itemgetter
import time, re, os, sys, gzip, pysam, copy, json, BedAnnoVar,BedAnnoAnno,Error
sys.dont_write_bytecode = True

class BedAnno(object):
	C3, C1, SO2Name, func2SO, Name2SO, Polar, C1toC3, AAnumber, AAcount, canonicalSS = dict(), dict(), dict(), dict(), dict(), dict(), dict(), dict(), None, dict()

	CURRENT_MT = 'NC_012920.1'

	C3 = {
		"AAA": "Lys", "AAC": "Asn", "AAG": "Lys", "AAT": "Asn",
		"ACA": "Thr", "ACC": "Thr", "ACG": "Thr", "ACT": "Thr",
		"AGA": "Arg", "AGC": "Ser", "AGG": "Arg", "AGT": "Ser",
		"ATA": "Ile", "ATC": "Ile", "ATG": "Met", "ATT": "Ile",
		"CAA": "Gln", "CAC": "His", "CAG": "Gln", "CAT": "His",
		"CCA": "Pro", "CCC": "Pro", "CCG": "Pro", "CCT": "Pro",
		"CGA": "Arg", "CGC": "Arg", "CGG": "Arg", "CGT": "Arg",
		"CTA": "Leu", "CTC": "Leu", "CTG": "Leu", "CTT": "Leu",
		"GAA": "Glu", "GAC": "Asp", "GAG": "Glu", "GAT": "Asp",
		"GCA": "Ala", "GCC": "Ala", "GCG": "Ala", "GCT": "Ala",
		"GGA": "Gly", "GGC": "Gly", "GGG": "Gly", "GGT": "Gly",
		"GTA": "Val", "GTC": "Val", "GTG": "Val", "GTT": "Val",
		"TAA": "*", "TAC": "Tyr", "TAG": "*", "TAT": "Tyr",
		"TCA": "Ser", "TCC": "Ser", "TCG": "Ser", "TCT": "Ser",
		"TGA": "*", "TGC": "Cys", "TGG": "Trp", "TGT": "Cys",
		"TTA": "Leu", "TTC": "Phe", "TTG": "Leu", "TTT": "Phe",

		"TCN": "Ser", "CCN": "Pro", "ACN": "Thr", "GTN": "Val",
		"CTN": "Leu", "GCN": "Ala", "CGN": "Arg", "GGN": "Gly",

		# inseq stop codon
		"UAA": "X", "UAG": "X",

		# selenocysteine
		"UGA": "Sec"
	}

	C1 = {
		"AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N",
		"ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
		"AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S",
		"ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I",
		"CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H",
		"CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
		"CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
		"CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
		"GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D",
		"GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
		"GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
		"GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
		"TAA": "*", "TAC": "Y", "TAG": "*", "TAT": "Y",
		"TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
		"TGA": "*", "TGC": "C", "TGG": "W", "TGT": "C",
		"TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F",

		"TCN": "S", "CCN": "P", "ACN": "T", "GTN": "V",
		"CTN": "L", "GCN": "A", "CGN": "R", "GGN": "G",

		"UAA": "X", "UAG": "X",

		"UGA": "U"
	}

	C1toC3 = dict()
	for threebase in sorted(C1.keys()):
		C1toC3[C1[threebase]] = C3[threebase]

	AAcount = len(C1toC3.keys())
	for k, v in zip(sorted(C1toC3.keys()) + ['?', '.'], range(1, AAcount + 3)):
		AAnumber[k] = v

	Polar = {
		'Ala': "NP", 'Arg': "P+", 'Asn': "P0", 'Asp': "P-",
		'Cys': "P0", 'Gln': "P0", 'Glu': "P-", 'Gly': "P0",
		'His': "P+", 'Ile': "NP", 'Leu': "NP", 'Lys': "P+",
		'Met': "NP", 'Phe': "NP", 'Pro': "NP", 'Sec': "NP",
		'Ser': "P0", 'Thr': "P0", 'Trp': "NP", 'Tyr': "P0",
		'Val': "NP",

		'X': '.', '*': '.'
	}

	for key in sorted(C1.keys()):
		Polar[C1[key]] = Polar[C3[key]]

		SO2Name = {

			# variant type
			"SO:0000159": 'del',
			"SO:1000032": 'delins',
			"SO:0001483": 'snv',
			"SO:0000667": 'ins',
			"ref": 'ref',
			"no-call": 'no-call',

			# Gene Parts
			"SO:0000316": 'CDS',
			"SO:0000204": 'five_prime_UTR',
			"SO:0000205": 'three_prime_UTR',
			"SO:0000655": 'ncRNA',
			"SO:0000191": 'interior_intron',
			"SO:0000448": 'three_prime_UTR_intron',
			"SO:0000447": 'five_prime_UTR_intron',
			"SO:0000163": 'five_prime_cis_splice_site',
			"SO:0000164": 'three_prime_cis_splice_site',
			"SO:0000167": 'promoter',
			"SO:0000605": 'intergenic_region',
			"span": 'span',

			# Function Parts
			"SO:0001819": 'synonymous_variant',
			"SO:0001583": 'missense_variant',
			"SO:0001587": 'stop_gained',
			"SO:0001578": 'stop_lost',
			"SO:0001822": 'inframe_deletion',
			"SO:0001821": 'inframe_insertion',
			"SO:0001589": 'frameshift_variant',
			"SO:0001582": 'initiator_codon_variant',
			"SO:0001893": 'transcript_ablation',
			"SO:0001567": 'stop_retained_variant',

			# for non-equal length of substitution
			"inframe_delins": 'inframe_delins',

			# for refSeq the same with call seq.
			"no-change": 'no-change',

			# for span and splice
			"unknown-likely-deleterious": 'unknown-likely-deleterious',

			# for ncRNA, utr, intron or intergenic
			"unknown": 'unknown',

			# for no-call variant
			"unknown-no-call": 'unknown-no-call',

			# the followings are replaced by 'unknown-likely-deleterious' in Voyager
			"SO:0001575": 'splice_donor_variant',
			"SO:0001574": 'splice_acceptor_variant',

			# newly added Functions in 1.01 for complementary of splicing variants
			"SO:0001572": 'exon_loss_variant',
			"SO:0001787": 'splice_donor_5th_base_variant',
			"SO:0001630": 'splice_region_variant',  # exon: 3bp + intron: 8bp
			"SO:0001995": 'extended_intronic_splice_region_variant',  # intron: 10bp

			# the followings are replaced by 'unknown' in Voyager
			"SO:0001623": '5_prime_UTR_variant',
			"SO:0001624": '3_prime_UTR_variant',
			"SO:0001619": 'nc_transcript_variant',
			"SO:0001627": 'intron_variant',

			"annotation-fail": 'annotation-fail',
			"abnormal-fs-site": 'abnormal-fs-site',
			"abnormal-intron": 'abnormal-intron',
			"abnormal-inseq-stop": "abnormal-inseq-stop",
		}

	Name2SO = {value: key for key, value in SO2Name.iteritems()}

	func2SO = {
		"abnormal-fs-site": 'abnormal-fs-site',
		"abnormal-inseq-stop": "abnormal-inseq-stop",
		"abnormal-intron": "abnormal-intron",
		"annotation-fail": 'annotation-fail',
		"cds-del": "SO:0001822",
		"cds-indel": "inframe_delins",
		"cds-ins": "SO:0001821",
		"cds-loss": "unknown-likely-deleterious",
		"coding-synon": "SO:0001819",
		"init-loss": "SO:0001582",
		"no-change": "no-change",
		"splice": "unknown-likely-deleterious",
		"splice-3": "unknown-likely-deleterious",
		"splice-5": "unknown-likely-deleterious",
		"stop-gain": "SO:0001587",
		"stop-loss": "SO:0001578",
		"stop-retained": "SO:0001567",
		"unknown-no-call": 'unknown-no-call',
		"utr-3": "unknown",
		"utr-5": "unknown",
		'altstart': "SO:0001582",
		'frameshift': "SO:0001589",
		'intron': "unknown",
		'knockout': "SO:0001893",
		'missense': "SO:0001583",
		'ncRNA': "unknown",
		'nonsense': "SO:0001587",
		'promoter': "unknown",
		'span': "unknown-likely-deleterious",
		'unknown': "unknown",

		# newly added in 1.01 for splicing variants complementary
		# may exists in alt_func keys, alt_funcSO, alt_funcSOname
		"exon-loss": "SO:0001572",
		"splice-5-5th": "SO:0001787",
		"splice-region": "SO:0001630",
		"splice-ext": "SO:0001995",
	}

	canonicalSS = {
		'D': "GT",
		'A': "AG",
	}

	GenePartsOrder = dict()
	gp = "CDS span five_prime_cis_splice_site three_prime_cis_splice_site ncRNA five_prime_UTR three_prime_UTR interior_intron five_prime_UTR_intron three_prime_UTR_intron abnormal-intron promoter annotation-fail intergenic_region".split()
	for i, j in zip(gp, range(1, 15)):
		GenePartsOrder[i] = j
	GenePartsOrder[""] = 15

	REF_BUILD = 'GRCh37'
	
	'''
	=head1 Methods

	=head2 new

	=over

	=item About : Creat a new annotation entry

	=item Usage :

		my $beda = BedAnno->new( db => "in.bed.gz", tr => 'in.trans.fas.gz', batch => 1 );

	=item Args    - (all database files should be tabix indexed)

	=over

	=item Essential Args:

	=over

	=item I<db> <in.bed.gz>

	=over

	=item annotation database file. 

	=back

	=item I<tr> <in.trans.fas.gz>

	=over

	=item transcript sequence fasta file

	=back

	=item See L</DATABASE FORMAT> for more infomation. 

	=back

	=item Optional Args :

	=over

	=item Common options :

	=over 

	=item I<cytoBand> [cytoBand.bed.gz]

	=over

	=item add cytoBand information

	=back

	=item I<rmsk> [rmsk.bed.gz]

	=over

	=item add repeat tag information

	=back

	=item I<gwas> [gwasCatalog_snp137.bed.gz]

	=over

	=item add gwas information depend on only position

	=back

	=item I<pfam> [pfam.tsv.gz]

	=over

	=item add pfam information

	=back

	=item I<prediction> [ensembl_prediction_db.tsv.gz]

	=over

	=item add sift, polyphen2 prediction and scores

	=back

	=item I<condel> [condel config path]

	=over

	=item compute condel scores and prediction based on sift and polyphen2 HumVar prediction

	=back

	=item I<cosmic> [Cosmic_v67_241013.bed.gz]

	=over

	=item add cosmic annotation information

	=back

	=item I<phyloP> [phyloP_scores.tsv.gz]

	=over

	=item add phyloP scores of all 3 datasets.

	=back

	=item I<dbSNP> [snp137.bed.gz]

	=over

	=item add rsID and dbSNP frequency information

	=back

	=item I<tgp> [tgp_phaseI_v3.bed.gz]

	=over

	=item add 1000 genomes allele frequency information

	=back

	=item I<cg54> [CG54.bed.gz]

	=over

	=item add 54 whole genomes allele frequency information from CompleteGenomics.

	=back

	=item I<wellderly> [wellderly.bed.gz]

	=over

	=item add wellderly's allele frequency information from CompleteGenomics.

	=back

	=item I<esp6500> [ESP6500.bed.gz]

	=over

	=item add ESP6500 allele frequency information from NHLBI

	=back

	=item I<exac> [ExAC.r0.2.sites.vep.vcf.gz]

	=over

	=item add ExAC allele frequency information from 61,486 unrelated individuals.

	=back

	=item I<customdb_XX> [custom db in the same format with esp6500's]

	=over

	=item add customdb XX frequency infomation, XX is the ID. (multiple db will require multiple customdb option)


	=back

	=item I<quiet>

	=over

	=item Suppress warning messege to output.

	=back

	=item I<batch> [boolean]

	=over

	=item use batch mode annotation, default in daemon mode as an annotation engine.

	=back

	=item I<genome> [ "refgenome.fa.rz" ]

	=over

	=item reference genome fasta, razipped and samtools faidxed for use.

	=back

	=item I<genes> [ "genes.list" | $rh_geneslist ]

	=over

	=item annotate transcripts for I<genes>. e.g. {"ABC" => 1, "DEF" => 1} or "genes.list" 

	=back

	=item I<trans> [ "trans.list" | $rh_translist ]

	=over

	=item annotate transcripts in I<trans>. e.g. {"NM_0012.1" => 1, "NM_0034.2" => 1} or "trans.list" 

	=back

	=item I<mmap> [boolean]

	=over

	=item allow annotating all other non "BEST" multiple-mapping records, boolen option, default not allowed.
		  e.g. NM_0123.1 have 3 mapping location or alternate splicing, default only the "BEST" one will be 
		  annotated. See L</DATABASE FORMAT> for the "BEST" definition.

	=back

	=back

	=item Batch mode options :

	=over

	=item I<region> [region_string]

	=over

	=item limit to only annotate transcript in I<region>. e.g. "chr20:1234567-1234568", prior to I<regbed>.

	=back
	 
	=item I<regbed> [BED format file]

	=over

	=item similar to I<region>, with allowing multiple regions. e.g. "in.region.bed". 

	=back

	=back

	=item Notes

	=over

	=item Batch mode is designed for annotation of a complete list of variants 
		  on same chromosome, read all information of the chr into memory, 
		  and annotate all variants together in the order of chr coordinates.
		  This mode can avoid frequent IO, brought by tabix searching, but need
		  huge memory cost.

	=back

	=back

	=back

	=item Returns

	=over

	=item Annotation Engine object entry, please see L</load_anno> for more information.

	=back

	=back

	=cut
	'''
	def __init__(self, kw):
		for k, v in kw.iteritems():
			setattr(self, k, v)

		if "db" not in kw or "tr" not in kw:
			raise Error.Error("Error: please at least give 'db' and 'tr' path.")

		debugOpt = 1 if "debug" in kw else 0
		t0 = None
		if debugOpt:
			t0 = [time.time()]

		if hasattr(self, "genes"):
			if type(self.genes) is not dict:
				try:
					GENE = open(self.genes, 'r')
				except IOError, e:
					raise IOError(self.genes + ":")
				genes = dict()
				for line in GENE:
					result, _ = re.subn(r"\s+", "", line)
					genes[result] = 1
				self.genes = genes
				GENE.close()
		if hasattr(self, "trans"):
			if type(self.trans) is not dict:
				try:
					TRAN = open(self.trans, 'r')
				except IOError, e:
					raise IOError(self.trans + ":")
				trans = dict()
				for line in TRAN:
					result, _ = re.subn(r"\s+", "", line)
					trans[result] = 1
				self.trans = trans
				TRAN.close()
			rclean_trans = dict()
			for k in self.trans.keys():
				result, _ = re.subn(r"\-\d+$", "", k)
				rclean_trans[result] = 1
			self.clean_trans = rclean_trans

		self.set_db(self.db)

		t1 = None
		if debugOpt:
			t1 = [time.time()]
			sys.stderr.write("{}{}{}".format("BedAnno->new [db load] ... ", str(t1 - t0), "\n"))
			t0 = t1

		self.set_tr(self.tr)
		if debugOpt:
			t1 = [time.time()]
			sys.stderr.write("{}{}{}".format("BedAnno->new [tr load] ... ", str(t1 - t0), "\n"))
			t0 = t1

		self.set_refbuild(BedAnno.REF_BUILD)

		for dbk in sorted(dir(self)):
			m = re.match(r'customdb_(\S+)', dbk)
			if m:
				dbID = m.group(1)
				self.set_customdb(eval("self." + dbk), dbID)  # set_customdb 函数没实现，稍后再转

		if debugOpt:
			t1 = [time.time()]
			sys.stderr.write("{}{}{}".format("BedAnno->new [db load] ... ", str(t1 - t0), "\n"))
			t0 = t1

	def set_db(self, db):
		if (not os.access(db, os.F_OK)) or (not os.access(db, os.R_OK)):
			raise Error.Error("Error: cannot read " + db + ".")

		setattr(self, 'db', db)
		if (not os.access("{}{}".format(db, ".tbi"), os.F_OK)) or (not os.access("{}{}".format(db, ".tbi"), os.R_OK)):
			raise Error.Error("Error: [" + db + "] index (.tbi) file not found, please build it first.")

		setattr(self, 'tidb', "data/test_db.bed.gz.tbi")

		if not hasattr(self, "batch"):
			return

		open_args = dict()
		if hasattr(self, "batch"):
			open_args["region"] = self.region
		if not hasattr(self, "region") and hasattr(self, "regbed"):
			open_args["regbed"] = self.regbed
		if hasattr(self, "genes"):
			open_args["genes"] = self.genes
		if hasattr(self, "trans") and hasattr(self, "clean_trans"):
			open_args["trans"] = self.trans
			open_args["clean_trans"] = self.clean_trans
		if hasattr(self, "mmap"):
			open_args["mmap"] = self.mmap

		self.annodb = self.load_anno(open_args)

	def throw(self, msg):
		print "".join(msg), "at", os.path.basename(__file__), "line ", sys._getframe().f_lineno
		raise Exception


	'''
	=head2 load_anno

		About   : load all needed annotation infomation into memory for multi-process annotation
		Usage   : my $ranndb = $beda->load_anno( region => "chr20:1234567-1234568", trans => \%trans );
		Args    : Using %args to override class's properties: region, regbed, genes, trans
				  if no args, use the the class's properties as default.
		Returns : a localized merged anno db, The returned annotation database is a hash ref.
			{
				$chr => [
					{
						sta   => $start, (0 based)
						sto   => $stop,  (1 based)
						annos => {
							$anno_string => $offset, ...
						}

						detail => {
							$tid => {
								gsym => $gsym,    (gene symbol)
								gid  => $gid,     (Entrez gene id)
								gpSO => $gpSO,    (GeneParts SO)
								blka => $blka,    (block attribute)
								exin => $exin,    (exon intron number)
								nsta => $nsta,    (n./r. left  of whole block)
								nsto => $nsto,    (n./r. right of whole block)
								csta => $csta,    (c.    left  of whole block)
								csto => $csto,    (c.    right of whole block)
								wlen => $wlen,    (length of whole block)
								pr   => $pr,      (primary tag)
								strd => $strd,    (strand)
								offset => $offset,(offset of current block to whole block)
								mismatch => $mismatch (non-equal length block descripter)
							}, ...
						}
					}, ... 
				], ...
			}
		  Note: when variation hit one of the annotation entry, the anno_string will be parsed.
		  and the "detail" tag will be added then.

	=cut

	'''
	def load_anno(self, args):
		query_region = list()
		if "region" in args:
			regions = re.split(r'\s+', args['region'])
			for reg in regions:
				if reg == "":
					continue
				m = re.search(r'^(\S+):(\-?\d+)\-(\d+)$', reg)
				if m:
					name, beg, end = m.group(1), m.group(2), m.group(3)
					if int(beg) <= 0:
						if hasattr(self, "debug"):
							self.warn(["Warning: region string should be 1 based [" + str(reg) + "], has been changed to 1 based"])
						beg = 1
					name = re.sub(r'^chr', "", name, flags=re.IGNORECASE)
					query_region.append([name, int(beg) - 1, end])
				else:
					self.throw(["Error: unavailable region string [" + str(reg) + "]."])
		elif ("regbed" in args) and os.path.exists(args['regbed']):
			BED = ""
			try:
				BED = open(args['regbed'], 'r')
			except IOError, e:
				self.throw(["Error: [" + args['regbed'] + "]"])
			for line in BED:
				line = line.strip('\r\n')
				beditm = re.split(r'\t', line)
				if 3 > len(beditm):
					self.throw(["Error: bed format error."])
				beditm[0] = re.sub(r'^chr', '', beditm[0], re.IGNORECASE)
				if re.match(r'^M', beditm[0], flags=re.IGNORECASE):
					beditm[0] = "MT"
				query_region.append(beditm[:3])
			BED.close()

		read_all_opt = 0
		all_querys = list()
		annodb_h = None
		if 0 == len(query_region):
			read_all_opt = 1
			try:
				annodb_h = gzip.open(self.db, 'rb')
			except:
				raise "Error: [" + self.db + "] " + "GunzipError\n"
		else:
			sorted_regions = sorted(query_region, key=itemgetter(0, 1, 2))
			cname, cbeg, cend = sorted_regions[0]
			tbx = pysam.TabixFile(self.db)
			for k in range(1, len(sorted_regions)):
				dname, dbeg, dend = sorted_regions[k]
				if dname != cname or dbeg > cend:
					try:
						query_ent = tbx.fetch(cname, cbeg, cend, parser=pysam.asTuple())
					except Exception, e:
						raise e
					all_querys.append(query_ent)
					cname, cbeg, cend = sorted_regions[k]
				else:
					cend = dend
			q = None
			try:
				q = tbx.fetch(cname, cbeg, int(cend),parser=pysam.asTuple())
			except Exception, e:
				raise e
			all_querys.append(q)
		# trans filter is always be ahead of genes
		if "mmap" in args:
			mmapTag = 1
		else:
			mmapTag = 0
		if "genes" in args:
			geneTag = 1
		else:
			geneTag = 0
		if "trans" in args:
			tranTag = 1
		else:
			tranTag = 0
		geneList = None
		tranList = None
		if geneTag:
			geneList = args['genes']
		if tranTag:
			tranList = args['trans']
		pureTran = dict()
		if tranTag:
			for k in tranList.keys():
				k = re.sub(r'\-\d+$', '', k)
				pureTran[k] = 1

		rannodb = dict()
		while True:
			tb_ent = None
			if read_all_opt:
				if type(annodb_h) is not file:
					break
				try:
					tb_ent = next(annodb_h)
				except:
					break
			else:
				while 0 < len(all_querys):
					region_read = None
					try:
						region_read = next(all_querys[0])
					except:
						del all_querys[0]
					else:
						tb_ent = region_read
						break
				if 0 == len(all_querys):
					break

			tb_ent = re.sub(r'\s+$', '', str(tb_ent))
			if len(re.split(r'\t', tb_ent)) < 4:
				self.throw("Error: db format unmatched")
			chr, start, stop, annostr = re.split(r'\t', tb_ent)[:4]
			if annostr == "":
				self.throw("Error: db format unmatched")

			annos = re.split(r'; ', annostr)

			ent = {"sta": start, "sto": stop}
			ent["annos"] = dict()
			for anno_ent in annos:
				cont = re.split(r'\|', anno_ent)
				# no multiple mapping transcripts if !$mmapTag
				if (re.search(r'\-\d+$', cont[0])) and (not mmapTag) and (not tranTag or cont[0] not in tranList):
					continue

				ori_cont = cont[0]
				cont[0] = re.sub(r'\-\d+$', '', cont[0])
				if (geneTag) and (cont[1] not in geneList) and (not tranTag or cont[0] not in pureTran):
					continue
				if (tranTag) and ((cont[0] not in pureTran) or (not mmapTag and ori_cont not in tranList)) and (
					not geneTag or cont[1] not in geneList):
					continue
				m = re.search(r'\|(\d+)$', anno_ent)

				ofst = None
				if m:
					ofst = m.group(1)
					anno_ent = re.sub(r'\|(\d+)$', '', anno_ent)
				else:
					self.throw("db format error: [" + anno_ent + "]")
				ent["annos"][anno_ent] = ofst
			if "annos" not in ent:
				continue
			if chr not in rannodb:
				rannodb[chr] = list()
			rannodb[chr].append(ent)

		if read_all_opt:
			annodb_h.close()

		if mmapTag or geneTag or tranTag:
			rannodb = BedAnno.region_merge(rannodb)

		return rannodb

	'''
	=head2 region_merge

		About   : merge consecutive same-entries regions
		Usage   : my $rannodb = region_merge($loaded_db);
		Args    : A hash ref of loaded_db.
		Returns : A hash ref of merged db.

	=cut
	'''
	@staticmethod
	def region_merge(radb):
		local_radb = copy.deepcopy(radb)

		for chr in local_radb.keys():
			annoents = copy.deepcopy(local_radb[chr])
			oricount = len(annoents)
			if oricount <= 1:
				continue
			# merge from bottom to top
			curid = oricount - 1
			while curid > 0:
				latter = annoents[curid]
				former = annoents[curid - 1]
				curid -= 1
				if latter['sta'] != former['sto']:
					continue
				former_ann = former['annos'].keys()
				formern = len(former_ann)
				if formern != len(latter['annos'].keys()):
					continue
				for ann in former_ann:
					if ann not in latter['annos']:
						break
				else:
					# splice the latter one and correct the former one
					former['sto'] = latter['sto']
					del annoents[curid + 1]
			if oricount > len(annoents):
				local_radb[chr] = copy.deepcopy(annoents)
		return local_radb

	def set_tr(self, tr):
		if (not os.access(tr, os.F_OK)) or (not os.access(tr, os.R_OK)):
			self.throw("Error: cannot read " + str(tr) + ".")
		self.tr = tr
		load_opts = dict()
		if hasattr(self, "genes"):
			load_opts['genes'] = self.genes
		if hasattr(self, "trans"):
			load_opts['trans'] = self.trans
		self.trInfodb = self.readtr(load_opts)

	'''
	=head2 readtr

		About   : Read transcript information, and even sequences if in batch mode.
		Usage   : my $rtrSeqs = $beda->readtr( genes => $rh_NewGenes, trans => $rh_NewTrans );
		Args    : Optional args "genes" and "trans" only accept hash ref values.
				  if no args specified, it will load information based on the
				  configuration of BedAnno entry.
		Returns : A hash ref of trSeqs:
				  {
					$tr_acc => {
						len      => $tr_len,
						gene     => $gene_sym,

						# optional tags:
						prot     => $prot_acc,
						plen     => $prot_len,
						csta     => $cds_start_on_trSeq, # 0 based
						csto     => $cds_end_on_trSeq,   # 1 based
						seq      => $tr_sequence,
						
				nfs      => { $fslead => $fsbase, ... },
				cfs      => { $cfslead => $fsbase, ... },
						X        => 1,                   # inseqStop
						U        => 1,                   # selenocysteine
						A        => 1,                   # polyATail
						altstart => {                    # altstart codons
							$startCodons1 => 1,
							$startCodons2 => 1,
							...
						},

						# the following two keys will be assigned
						# to mRNA when needed
						cseq     => $codonSequence,
						pseq     => $proteinSequence,
					},
					...
				  }

	=cut
	'''
	def readtr(self, opts):
		if "genes" in opts and type(opts["genes"]) is not dict:
			self.throw("Options arg 'genes' only accept hash ref as value.")
		if "trans" in opts and type(opts["trans"]) is not dict:
			self.throw("Options arg 'trans' only accept hash ref as value.")

		fas_h = None
		if re.search(r'\.gz$', self.tr):
			try:
				fas_h = gzip.open(self.tr, 'rb')
			except Exception, e:
				self.throw("Error: [" + str(self.tr) + "] GunzipError+\n")
		else:
			try:
				fas_h = open(self.tr)
			except Exception, e:
				self.throw("Error: [" + str(self.tr) + "] \n")
		seqs = dict()
		for line in BedAnno.fasta_line(fas_h, ">"):
			line = re.sub(r'[\s>]+$', '', line)
			if re.search(r'^\s*$', line):
				continue
			hd=None
			m=re.search(r'^(\S+[^\n]*)\n',line)
			if m:
				hd =m.group(1)
				line=re.sub(r'^(\S+[^\n]*)\n','',line)
			#hd = BedAnno.re_sub_b_ref(line, "^(\S+[^\n]*)\n")
			if hd is None:
				self.throw("Error: trSeq parse error!")
			headers = re.split(r'\s+', hd, 10)
			if len(headers) < 7:
				self.throw("Error: trSeq header parse error!")

			line = re.sub(r'\s+', '', line)
			if "trans" not in opts and "genes" not in opts:
				if not ((not hasattr(self, "clean_trans") and not hasattr(self, "genes")) or (
					hasattr(self, "clean_trans") and headers[0] in self.clean_trans) or (
					hasattr(self, "genes") and headers[2] in self.genes)):
					continue
			else:
				if not (("trans" in opts and headers[0] in opts['trans']) or (
						"genes" in opts and headers[2] in opts['genes'])):
					continue
			seqs[headers[0]] = dict()
			seqs[headers[0]]['seq'] = line
			seqs[headers[0]]['len'] = headers[1]  # tx length
			seqs[headers[0]]['gene'] = headers[2]  # gene symbol
			if headers[3] != ".":  # prot acc.ver
				seqs[headers[0]]['prot'] = headers[3]
			if headers[4] != ".":  # prot length
				seqs[headers[0]]['plen'] = headers[4]
			if re.search(r'selenocysteine', headers[6]):  # selenocysteine
				seqs[headers[0]]['U'] = 1
			if re.search(r'inseqStop', headers[6]):  # inseqStop
				seqs[headers[0]]['X'] = 1
			if re.search(r'polyATail', headers[6]):  # polyATail
				seqs[headers[0]]['A'] = 1
			if headers[5] != ".":  # cds start, cds end in tx
				seqs[headers[0]]['csta'] = re.split(r',', headers[5])[0]
				seqs[headers[0]]['csto'] = re.split(r',', headers[5])[1]
			fss = dict()

			if re.search(r'altstart', headers[6]):
				if 'altstart' not in seqs[headers[0]]:
					seqs[headers[0]]['altstart']=dict()
				for i in re.split(r';', headers[7]):
					seqs[headers[0]]['altstart'].update({i: 1})
			elif len(headers) >= 8:  # frameshift at 8th item
				fss = re.split(r"\|", headers[7])
			if len(headers) >= 9:  # frameshift at 9th item
				fss = re.split(r"\|", headers[8])
			if len(fss) > 0:
				if headers[5] == ".":
					self.throw("Error: [" + str(headers[0]) + "] no cds but with fs!")
				seqs[headers[0]]['nfs'] = dict(map(lambda x: (x.split(";")[0], x.split(";")[1]), fss))
				seqs[headers[0]]['cfs'] = dict(
					map(lambda x: (int(x.split(";")[0]) - int(seqs[headers[0]]['csta']), x.split(";")[1]), fss))
		fas_h.close()
		return seqs

	@staticmethod
	def fasta_line(fh, sep):
		s_out = ""
		for s_v in fh:
			if sep in s_v:
				yield s_out + s_v[:(s_v.index(sep) + 1)]
				s_out = s_v[(s_v.index(sep) + 1):]
			else:
				s_out = s_out + s_v
		yield s_out

	@staticmethod
	def re_sub_b_ref(line, reg):
		m = re.search(reg, line)
		hd = None
		if m:
			hd = m.group(1)
		return hd

	'''
	=head2 set/get methods for properties

		List of Properties:
							get          set
		db                  o            o
		tr                  o            o
		refbuild            o            o
		tidb                o            x
		annodb              o            x
		trInfodb            o            x
		cytoBand            o            o
		cytoBand_h          o            x
		rmsk                o            o
		rmsk_h              o            x
		gwas                o            o
		gwas_h              o            x
		pfam                o            o
		pfam_h              o            x
		prediction          o            o
		prediction_h        o            x
		phyloP              o            o
		phyloP_h            o            x
		cosmic              o            o
		cosmic_h            o            x
		dbSNP               o            o
		dbSNP_h             o            x
		tgp                 o            o
		tgp_h               o            x
		esp6500             o            o
		esp6500_h           o            x
		exac                o            o
		exac_h              o            x
		cg54                o            o
		cg54_h              o            x
		wellderly           o            o
		wellderly_h         o            x

		e.g.    : $beda->set_refbuild($refbuild);
				  my $refbuild = $beda->get_refbuild();

	=cut
	'''
	def set_refbuild(self, custom_build_info):
		self.refbuild = custom_build_info  # 稍后再转

	'''
	=head2 varanno

		About   : implicitly create a new BedAnno::Anno entry, and 
				  assign all the needed annotation for to it.
		Usage   : ($rAnnoRst, $AEIndex) = $beda->varanno($var, $AEIndex);
		Args    : The BedAnno entry, BedAnno::Var entry and current dbidx, 
				  current dbidx should be always 0, if used for non-batch mode.
		Returns : A BedAnno::Anno entry and current dbidx for nex query in batch.
				{
					var => {

						# the first part is from var parsing result.
						# please see "BedAnno::Var new()".
						# import all the keys from original parsed var entry
						# and add the following keys by this method.

						varName => $var_mutation_name,

						# information
						varTypeSO => $varTypeSO,
						gHGVS     => $gHGVS,
						refbuild  => $referenceBuild,

						# Here's some optional parts which may be generated
						# when extra resource is available:

						cytoBand  => $cytoBand,
						reptag    => $repeatTag,
						gwas      => $ref_gwas_ret,

						# For single position for now
						phyloPpm    => $PhyloPscorePlacentalMammals,
						phyloPpr    => $PhyloPscorePrimates,
						phyloPve    => $PhyloPscoreVetebrates,

						reptag	=> $RepeatMaskerTag,

						dbsnp => {
							$rsID => {
								AN => $dbsnp_total_allele_count,
								AF => $dbsnp_alt_allele_frequency,  # though ref
							},
							...
						},

						cosmic => $ref_cosmic_return,

						tgp => {
							AN => $tgp_total_allele_count,
							AF => $tgp_alt_allele_frequency,
						},

						cg54 => {
							AN => $cg54_total_allele_count,
							AF => $cg54_alt_allele_frequency,
						},

						wellderly => {
							AN => $wellderly_total_allele_count,
							AF => $wellderly_alt_allele_frequency,
						},

						esp6500 => {
							AN => $esp6500_total_allele_count,
							AF => $esp6500_alt_allele_frequency,
						},

						cusdb_XX => {
						AN => $custom_db_allele_count,
						AF => $custom_db_allele_frequency,
						},
						...

					},
					trInfo => {
						$tid => {
							trVarName     => $transcriptVariantName,
							geneId        => $Entrez_Gene_ID,
							geneSym       => $Gene_Symbol,
							prot          => $Protein_Acc_Ver,
							strd          => $strand,
							rnaBegin      => $Begin_in_RNA_transcript,
							rnaEnd        => $End_in_RNA_transcript,
							cdsBegin      => $Begin_in_CDS,            # cDot format
							cdsEnd        => $End_in_CDS,              # cDot format
							protBegin     => $Begin_in_Protein,
							protEnd       => $End_in_Protein,
							c             => $cHGVS,
							p             => $pHGVS,
							p3	          => $threeletter_pHGVS,
							cc            => $codon_change,
							polar         => $polar_change,
							r             => $imp_funcRegion,
							r_Begin       => $imp_beginfuncRegion,
							r_End         => $imp_endfuncRegion,
							func          => $imp_funcCode,
							exin          => $exIntr_number,
							ei_Begin      => $imp_Begin_exIntr_number,
							ei_End        => $imp_End_exIntr_number,
							genepart      => $GenePart,
							genepartSO    => $GenePartSO,
							componentIndex => $componentIndex,
							exonIndex     => $exonIndex,               # '.' for N/A
							intronIndex   => $intronIndex,             # '.' for N/A
							funcSOname    => $FunctionImpact,
							funcSO        => $FunctionImpactSO,
							trAlt         => $alt_string_on_transcript,
							trRef         => $ref_string_on_transcript,
							prAlt         => $protein_alt_sequence,
							prRef         => $protein_ref_sequence,
							primaryTag    => $refstandard_primary_or_not,	# Y/N
							preStart => {    # the position before the start of var
								nDot => $rna_hgvs_pos,
								cDot => $cds_hgvs_pos,
								r    => $func_region,
								exin => $exon_intron_number,
							},
							postEnd => {     # the position after the end of var
								nDot => $rna_hgvs_pos,
								cDot => $cds_hgvs_pos,
								r    => $func_region,
								exin => $exon_intron_number,
							},
							trRefComp => {

								# some trRef components
							},

							# for some of splice variants, there may exists 
							# the following function information
							alt_func        => $alternative_func_code,
							alt_funcSO      => $alternative_variant_SO_id,
							alt_funcSOname  => $alt_variant_SO_name,

							# The following parts will be exists if extra resource
							# is available.
							pfamId      => $PFAM_ID,
							pfamName    => $PFAM_NAME,
							siftPred    => $SIFTpred,
							siftScore   => $SIFTscore,
							pp2divPred  => $Polyphen2HumDivPred,
							pp2divScore => $Polyphen2HumDivScore,
							pp2varPred  => $Polyphen2HumVarPred,
							pp2varScore => $Polyphen2HumVarScore,
							condelPred  => $Condelpred,
							condelScore => $Condelscore,

						},
						...
					}
				}

	=cut
	'''
	def varanno(self, var, AEIndex=None):
		if AEIndex is None:
			AEIndex = 0

		debugOpt = None
		if hasattr(self, "debug"):
			debugOpt = 1
		else:
			debugOpt = 0

		if debugOpt:
			t0 = [time.time()]

		if hasattr(self, "cytoBand"):
			var.cytoBand = self.cytoBand_h.getCB(var.chr, var.pos, var.end)  # 注意调用函数的实现

		if hasattr(self, "rmsk"):
			var.reptag = self.rmsk_h.getRepTag(var.chr, var.pos, var.end)  # 注意调用函数的实现

		if hasattr(self, "gwas"):
			var.gwas = self.gwas_h.getGWAS(var.chr, var.pos, var.end)  # 注意调用函数的实现

		if hasattr(self, "phyloP"):
			if var.sm == 1:
				var.phyloPpm, var.phyloPpr, var.phyloPve = self.phyloP_h.getPhyloP46wayScore(var.chr, (
				int(var.pos) + 1))  # 注意调用函数的实现

		if hasattr(self, "dbSNP"):
			if hasattr(var, "sep_snvs"):
				new_sqls = dict()
				cur_start = var.sep_snvs[0]
				for i in range(0, len(var.sep_snvs)):
					if i == (len(var.sep_snvs) - 1) or var.sep_snvs[i + 1] - cur_start > 1:
						new_ref = var.ref[(new_ref - int(var.pos) - 1):(int(var.sep_snvs[i]) - var.pos - 1)]
						new_alt = var.alt[(cur_start - int(var.pos) - 1):(int(var.sep_snvs[i]) - var.pos - 1)]
						new_sqls.append([cur_start, var.sep_snvs[i], new_ref, new_alt])
						if i < len(var.sep_snvs) - 1:
							cur_start = var.sep_snvs[i + 1]
				for rSE in new_sqls:
					rOneSql = self.dbSNP_h.getRS(var.chr, rSE)  # 注意调用函数的实现
					for k in sorted(rOneSql.keys()):
						var.dbsnp[k] = rOneSql[k]
			else:
				var.dbsnp = self.dbSNP_h.getRS(var.chr, var.pos, var.end, var.ref, var.alt)  # 注意调用函数的实现

		if hasattr(self, "tgp"):
			var.tgp = self.tgp_h.getAF(var.chr, var.pos, var.end, var.ref, var.alt)  # 注意调用函数的实现

		if hasattr(self, "esp6500"):
			var.esp6500 = self.esp6500_h.getAF(var.chr, var.pos, var.end, var.ref, var.alt)  # 注意调用函数的实现

		if hasattr(self, "exac"):
			var.exac = self.exac_h.getAF(var.chr, var.pos, var.end, var.ref, var.alt)  # 注意调用函数的实现

		for dbhk in sorted(dir(self)):
			m = re.match(r'cusdb_(\S+)_h', dbhk)
			if m and hasattr(self, dbhk):
				dbID = m.group(1)
				setattr(var, "cusdb_" + dbID, self.dbhk.getAF(var.chr, var.pos, var.end, var.ref, var.alt))
				# eval("var."+"cusdb_"+dbID)=self.dbhk.getAF(var.chr,var.pos,var.end,var.ref,var.alt)#注意调用函数的实现

		if hasattr(self, "cg54"):
			var.cg54 = self.cg54_h.getAF(var.chr, var.pos, var.end, var.ref, var.alt)  # 注意调用函数的实现

		if hasattr(self, "wellderly"):
			var.wellderly = self.wellderly_h.getAF(var.chr, var.pos, var.end, var.ref, var.alt)  # 注意调用函数的实现

		if hasattr(self, "cosmic"):
			var.cosmic = self.cosmic_h.getCOSMIC(var.chr, var.pos, var.end, var.ref, var.alt)  # 注意调用函数的实现

		t1 = None
		if debugOpt:
			t1 = [time.time()]
			sys.stderr.write("BedAnno->anno [var extradb sql] ... " + str(t1 - t0) + "\n")
			t0 = t1
		var.varTypeSO = BedAnno.Name2SO[var.guess]
		var.refbuild = self.refbuild
		var.get_gHGVS()
		if re.search(r'\]$', var.gHGVS):
			var.standard_gHGVS = BedAnno.gen_standard_gphgvs(var.gHGVS)
			var.alt_gHGVS = BedAnno.gen_alt_ghgvs(var.gHGVS)

		# Due to the bed format database,
		# add flank left and right 1bp to query the database
		# to involve all the mismatches

		var_region = None
		tmp = None
		if var.pos > 0:
			tmp = int(var.pos) - 1
		else:
			tmp = 0
		var_region = var.chr + ":" + str(tmp) + '-' + str(int(var.pos) + int(var.reflen) + 1)
		localdb = None
		if not hasattr(self, "batch"):  # daemon mode
			locOpts = {"region": var_region}
			if hasattr(self, "trans"):
				locOpts['trans'] = self.trans
			if hasattr(self, "genes"):
				locOpts['genes'] = self.genes
			locLoad = self.load_anno(locOpts)
			if var.chr in locLoad:
				localdb = locLoad[var.chr]
			else:
				localdb = []
		else:
			if var.chr in self.annodb:
				if debugOpt:
					self.warn("Warning: incomplete annotation database in batch mode, [" + str(var.chr) + "]")
				localdb = []
			else:
				localdb = self.annodb[var.chr]

		if debugOpt:
			t1 = [time.time()]
			sys.stderr.write("BedAnno->varanno [localdb load] ... " + str(t1 - t0) + "\n")
			t0 = t1

		annoEnt = BedAnnoAnno.BedAnnoAnno(var)

		if debugOpt:
			t1 = [time.time()]
			sys.stderr.write("BedAnno->varanno [new BedAnno::Anno] ... " + str(t1 - t0) + "\n")
			t0 = t1

		AEIndex = annoEnt.getTrPosition(localdb, AEIndex)
		if debugOpt:
			t1 = [time.time()]
			sys.stderr.write("BedAnno->varanno [getTrPosition] ... " + str(t1 - t0) + "\n")
			t0 = t1

		annoEnt = self.getTrChange(annoEnt)
		if debugOpt:
			t1 = [time.time()]
			sys.stderr.write("BedAnno->varanno [getTrChange] ... " + str(t1 - t0) + "\n")
			t0 = t1

		annoEnt = self.finaliseAnno(annoEnt)
		if debugOpt:
			t1 = [time.time()]
			sys.stderr.write("BedAnno->varanno [finaliseAnno] ... " + str(t1 - t0) + "\n")
			t0 = t1

		return (annoEnt, AEIndex)

	'''
	=head2 getTrChange

		About   : Calculate the transcript changes, based on TrPostition
		Usage   : $beda->getTrChange($annoEnt);
		Returns : assign the following tags in annoEnt
					trRef, prot, c, p, cc, polar, func
					prRef, prAlt

	=cut
	'''
	def getTrChange(self, annoEnt):
		# do nothing if not hit any transcript
		if not hasattr(annoEnt,"trInfo") or (0 == len(annoEnt.trInfo.keys())):
			return annoEnt

		getseq_cache = dict()
		for tid in sorted(annoEnt.trInfo.keys()):
			trannoEnt = annoEnt.trInfo[tid]
			unify_p, unify_r, unify_a, unify_rl, unify_al = annoEnt.var.getUnifiedVar(trannoEnt['strd'])

			# 1. check if annotation-fail
			if str(trannoEnt['genepartSO']) == 'annotation-fail' or str(trannoEnt['rnaBegin']) == '?' or str(trannoEnt['rnaEnd']) == '?':
				trannoEnt['func'] = 'annotation-fail'
				continue

			# the database hash of tr don't hash the annotation failed transcript
			qtid = tid
			qtid = re.sub(r'\-\d+$', '', str(qtid)) # trim the multiple mapping indicator in trAcc
			if qtid not in self.trInfodb or str(self.trInfodb[qtid]['seq']) == "":
				self.throw("".format("Error: your fasta database file may not complete. [", qtid, "]"))
			trdbEnt = self.trInfodb[qtid]
			trSeq = trdbEnt['seq']
			if "cseq" not in trdbEnt and "csta" in trdbEnt:
				trdbEnt['cseq'] = BedAnno._getCodingSeq(trdbEnt)
			cdsOpt = None
			if "prot" in trdbEnt and trdbEnt['prot'] != "." and trannoEnt['cdsBegin'] != '':
				cdsOpt = 1
			else:
				cdsOpt = 0
			strd = None
			if trannoEnt['strd'] == '+':
				strd = 1
			else:
				strd = 0

			# fix the trRefComp when ended in 3'downstream

			# debug
			# print STDERR Data::Dumper->Dump(
			#    [$tid, $trannoEnt, $unify_r, $trSeq, $strd],
			#    ["tid", "trannoEnt", "unify_r", "trSeq", "strd"] );

			trRef = BedAnno.getTrRef(trannoEnt, unify_r, trSeq, strd)
			trannoEnt['trRef'] = trRef

			trAlt = trannoEnt['trAlt']

			if "prot" in trdbEnt and trdbEnt['prot'] != ".":
				trannoEnt['prot'] = trdbEnt['prot']
			f = 'c.' if cdsOpt else 'n.'
			chgvs_5 = trannoEnt['cdsBegin'] if cdsOpt else trannoEnt['rnaBegin']
			chgvs_3 = trannoEnt['cdsEnd'] if cdsOpt else trannoEnt['rnaEnd']
			if cdsOpt:
				trannoEnt['protBegin'], trannoEnt['protEnd'] = BedAnno._getCoveredProd(trdbEnt, chgvs_5, chgvs_3)
			# [ aaPos, codon, aa, polar, frame, [framealt] ]
			rcInfo5, rcInfo3 = BedAnno._getPairedCodon(trdbEnt, trannoEnt['rnaBegin'], trannoEnt['rnaEnd'])

			# we just use position comparison to show transcript ref
			# stat, instead of sm in var, because there may exists
			# indel difference bewteen refgenome and refSeq
			cmpPos = self.cmpPos(chgvs_5, chgvs_3)

			# 2. check if no-call
			if str(trAlt) == '?':
				trannoEnt['func'] = 'unknown-no-call'
				if int(cmpPos) == 0:  # 1 bp
					trannoEnt['c'] ="{}{}{}{}".format( f , chgvs_5 , trRef , '>?')
				elif int(cmpPos) == 1:  # multiple bp
					tmp = None
					m = re.search(r'^[ACGTN]+$', str(trRef))
					if m:
						tmp = trRef
					else:
						tmp = ""
					trannoEnt['c'] = "{}{}{}{}{}{}{}".format(f, chgvs_5, '_', chgvs_3, 'del', tmp, 'ins?')
				else:  # ins : reverse the positions
					trannoEnt['c'] = "{}{}{}{}{}".format(f, chgvs_3, '_', chgvs_5, 'ins?')

				if int(rcInfo5[0]) == int(rcInfo3[0]) and rcInfo5[0] > 0:
					# single codon
					aaOut = rcInfo5[2]
					trannoEnt['p'] = "{}{}{}{}".format('p.', aaOut, rcInfo5[0], '?')
					trannoEnt['cc'] = "{}{}".format(rcInfo5[1], '=>?')
					trannoEnt['polar'] = "{}{}".format(rcInfo5[3], '=>?')
				elif int(rcInfo3[0]) > 0 and int(rcInfo5[0]) > 0:  # multiple codon

					aa5 = rcInfo5[2] + rcInfo5[0]
					aa3 = rcInfo3[2] + rcInfo3[0]

					diopt = None
					if int(rcInfo5[0]) < int(rcInfo3[0]):
						diopt = 1
					else:
						diopt = 0

					if diopt:
						trannoEnt['p'] = "{}{}{}{}{}{}".format('p.', aa5, '_', aa3, 'del', 'ins?')
					else:
						trannoEnt['p'] = "{}{}{}{}{}{}".format('p.', aa3, '_', aa5, '', 'ins?')
				continue
			if re.search(r'=', str(trRef)):  # reference sequence too long
				if trannoEnt['r_Begin'] != trannoEnt['r_End']:
					trannoEnt['func'] = 'span'
				else:
					trannoEnt['func'] = 'unknown'
				trannoEnt['c'] = "{}{}{}{}{}".format(f, chgvs_5, '_', chgvs_3, 'del')
				if trAlt != "":
					trannoEnt['c'] = "{}{}{}".format(trannoEnt['c'], 'ins', trAlt)
				continue
			if str(trRef) == str(trAlt):
				trannoEnt['func'] = 'no-change'
				trannoEnt['c'] = f
				if int(cmpPos) == 0:
					trannoEnt['c'] = "{}{}{}{}".format(trannoEnt['c'], chgvs_5, trRef, '=')
				elif int(cmpPos) > 0:
					trannoEnt['c'] = "{}{}{}{}{}".format(trannoEnt['c'], chgvs_5, "_", chgvs_3, "=")
				else:
					# hgvs 2.1511 don't refer to this kind of no-change, which
					# come from del/ins mismatch from refgenome and refseq
					# with 0 length refSeq.
					trannoEnt['c'] = "{}{}".format(trannoEnt['c'], "=")
				continue

			trBegin, trEnd, real_var, rUnified = self.trWalker(tid, trannoEnt)
			real_p, real_r, real_a, real_rl, real_al = rUnified

			if cdsOpt:
				chgvs_5 = BedAnno._cPosMark(trBegin, trdbEnt['csta'], trdbEnt['csto'], trdbEnt['len'])
			else:
				chgvs_5 = trBegin
			if cdsOpt:
				chgvs_3 = BedAnno._cPosMark(trEnd, trdbEnt['csta'], trdbEnt['csto'], trdbEnt['len'])
			else:
				chgvs_3 = trEnd
			cmpPos = self.cmpPos(chgvs_5, chgvs_3)

			# check if a rep get across the utr/cds edge and can be curated
			# to a variant in cds region, add a new key uncurated_cHGVS
			# to restore the uncurated version of cHGVS
			if str(real_var.imp) == 'rep' and cdsOpt and (bool(re.search(r'^\d+$', str(chgvs_5))) != bool(re.search(r'^\d+$', str(chgvs_3)))) and bool(re.search(r'^\-\d+$', str(chgvs_5))) != bool(re.search(r'^\*\d+$', str(chgvs_3))):
				lesser_len = real_rl if int(real_rl) < int(real_al) else real_al
				m1 = re.search(r'^\-(\d+)$', str(chgvs_5))
				m2 = re.search(r'^(\d+)$', str(chgvs_5))
				if (m1 and int(m1.group(1)) <= int(lesser_len)) or (m2 and (int(trdbEnt['csto']) - int(m2.group(1))) <= int(lesser_len)):
					# 5utr and first cds
					mm=None
					if m1:
						mm=m1
					else:
						mm=m2
					offset = mm.group(1)
					opt_53 = 1 if re.search(r'^\-', str(chgvs_5)) else 0
					if not opt_53:
						offset = int(trdbEnt['csto']) - int(offset)
					change_cn = 0
					while int(offset) > int(change_cn) * int(real_var.replen):
						change_cn += 1
					# restore uncurated version of cHGVS
					trannoEnt['uncurated_cHGVS'] = "{}{}{}{}{}".format(f, chgvs_5, '_', chgvs_3, 'del')
					if re.search(r'^[ACGTN]+$', str(real_r)) and len(real_r) < 50:
						trannoEnt['uncurated_cHGVS'] = "{}{}".format(trannoEnt['uncurated_cHGVS'], real_r)
					else:
						trannoEnt['uncurated_cHGVS'] = "{}{}".format(trannoEnt['uncurated_cHGVS'], "")

					if int(real_al) > 0:
						trannoEnt['uncurated_cHGVS'] = "{}{}{}".format(trannoEnt['uncurated_cHGVS'], 'ins' , real_a)

					trPosChanged = int(change_cn) * int(real_var.replen)
					trBegin += trPosChanged
					real_var = BedAnnoVar.BedAnnoVar(tid, 0, real_rl - trPosChanged, real_r[trPosChanged:], real_a[trPosChanged:])

					real_p, real_r, real_a, real_rl, real_al = real_var.getUnifiedVar('+')
					if cdsOpt:
						chgvs_5 = BedAnno._cPosMark(trBegin, trdbEnt['csta'], trdbEnt['csto'], trdbEnt['len'])
					else:
						chgvs_5 = trBegin
					if cdsOpt:
						chgvs_3 = BedAnno._cPosMark(trEnd, trdbEnt['csta'], trdbEnt['csto'], trdbEnt['len'])
					else:
						chgvs_3 = trEnd
					cmpPos = self.cmpPos(chgvs_5, chgvs_3)

			# debug
			#        print STDERR Data::Dumper->Dump(
			#            [ $real_var, $trBegin,  $trEnd,  $chgvs_5,  $chgvs_3,  $cmpPos ],
			#            [ "real_var", "trBegin", "trEnd", "chgvs_5", "chgvs_3", "cmpPos" ]
			#        );

			# * check if a repeat case, and use cerntain chgvs string
			# * assign cHGVS string
			# not ( bool(not re.serach(r'^\d+$',chgvs_5))!=bool(not re.serach(r'^\d+$',chgvs_3)))      not get across the edge of cds/non-cds region
			# not(bool(not re.serach(r'\d+[\+\-][ud]?\d+',chgvs_5))!=bool(not re.serach(r'\d+[\+\-][ud]?\d+',chgvs_3)))     not get across the edge of intron/exon region or promoter / 3'd
			if real_var.imp == 'rep' and not (
				bool(not re.search(r'^\d+$', str(chgvs_5))) != bool(not re.search(r'^\d+$', str(chgvs_3)))) and not (
				bool(not re.search(r'\d+[\+\-][ud]?\d+', str(chgvs_5))) != bool(
					not re.search(r'\d+[\+\-][ud]?\d+', str(chgvs_3)))):
				# we currently assume repeat variant won't get across more than
				# two region, so ommit other case test
				trRep = real_var.rep
				if int(real_var.ref_cn) == 1 and int(real_var.alt_cn) == 2:
					if int(real_var.replen) == 1:
						trannoEnt['c'] = "{}{}{}{}".format(f, chgvs_5, 'dup', trRep)
					else:
						trannoEnt['c'] = "{}{}{}{}{}{}".format(f, chgvs_5, '_', chgvs_3, 'dup', trRep)
				else:
					trannoEnt['c'] = "{}{}{}{}{}{}{}{}".format(f, chgvs_5, trRep, '[', real_var.ref_cn, '>',
															   real_var.alt_cn, ']')

				# 0.4.9
				# add a new key: alt_cHGVS, using for query database which
				# do not have current standard cHGVS string, especially
				# for repeat case
				#
				# 0.5.0
				# generate standard_cHGVS for this kind of variants
				# prior to alt_cHGVS, change alt_cHGVS to only ins/del
				# format notation
				#
				# give this opt to indicate whether to give standard_cHGVS
				# a duplication format notation
				dupopt = None
				if int(real_var.ref_cn) > 1 and int(real_var.alt_cn) - int(real_var.ref_cn) == 1:
					dupopt = 1
				else:
					dupopt = 0
				if int(real_var.ref_cn) > int(real_var.alt_cn) or dupopt:
					# deletion and duplication
					if not dupopt:
						trannoEnt['standard_cHGVS'] = "{}{}{}{}{}{}".format(f, chgvs_5, trRep, '[', real_var.alt_cn,
																			']')
					changed_cn = abs(real_var.ref_cn - real_var.alt_cn)
					changed_type = 'dup' if dupopt else 'del'
					changed_cont = trRep * changed_cn
					changed_len = len(changed_cont)
					if int(changed_cn) == 1 and int(real_var.replen) == 1:
						single_change = "{}{}{}{}".format(f, chgvs_3, changed_type, trRep)
						if dupopt:
							trannoEnt['standard_cHGVS'] = single_change
						else:
							trannoEnt['alt_cHGVS'] = single_change
					else:
						renew_offset = real_var.replen * (
						real_var.ref_cn - 1) if dupopt else real_var.replen * real_var.alt_cn
						renew_cHGVS = None
						m1 = re.search(r'^([\*\+]?)(\-?\d+)$', str(chgvs_5))
						m2 = re.search(r'^(\-?[^\-\+]+)(\+?d?)(\-?u?\d+)$', str(chgvs_5))
						if re.search(r'^\d+$', str(chgvs_5)):  # cds / ncRNA exon
							renew_cHGVS = "{}{}{}{}{}{}".format(f, chgvs_5 + renew_offset, '_', chgvs_3, changed_type,
																changed_cont)
						elif m1:
							# utr region / ncRNA promoter 3'd region
							sig = m1.group(1)
							start_pos = m1.group(2)
							if sig == '':  # 5 utr
								renew_cHGVS = "{}{}{}{}{}{}".format(f, start_pos + renew_offset, '_', chgvs_3,
																	changed_type, changed_cont)
							else:  # 3 utr
								renew_cHGVS = "{}{}{}{}{}{}{}".format(f, sig, start_pos + renew_offset, '_', chgvs_3,
																	  changed_type, changed_cont)
						elif m2:
							# intron or promoter or 3'downstream
							anchor5 = m2.group(1)
							sig5 = m2.group(2)
							offset5 = m2.group(3)
							u5_opt = 0
							if re.search(r'u', offset5):
								u5_opt = 1
								offset5 = re.sub(r'u', '', offset5)
							m = re.search(r'^(\-?[^\-\+]+)(\+?d?)(\-?u?\d+)$', str(chgvs_3))
							if m:
								anchor3 = m.group(1)
								sig3 = m.group(2)
								offset3 = m.group(3)
								u3_opt = 0
								if re.search(r'u', str(offset3)):
									u3_opt = 1
									offset3 = re.sub(r'u', '', str(offset3))
								if anchor3 == anchor5 and sig5 == sig3:
									anc_5 = int(offset5) + int(renew_offset)
									if u5_opt:
										anc_5 = re.sub(r'^-', '-u', anc_5)
									renew_cHGVS = "{}{}{}{}{}{}{}{}".format(f, anchor5, sig5, anc_5, '_', chgvs_3,
																			changed_type, changed_cont)
								elif str(anchor3) != str(anchor5) and str(sig5) == '+' and str(sig3) == '':
									# assume same intron
									half_offset = (int(real_rl) - (int(offset5) + int(offset3))) / 2
									if int(renew_offset) < int(half_offset):
										renew_cHGVS = "{}{}{}{}{}{}{}{}".format(f, anchor5, sig5,
																				int(offset5) + int(renew_offset), '_', chgvs_3,
																				changed_type, changed_cont)
									else:
										renew_cHGVS = "{}{}{}{}{}{}{}{}".format(f, anchor3, sig3,
																				int(offset3) - int(changed_len) + 1, '_', chgvs_3,
																				changed_type, changed_cont)
								else:
									if not hasattr(self, "quiet"):
										self.warn("{}{}{}{}".format(
											"Warning: repeat get across more than 3 different region: ", tid, ":",
											trannoEnt['c']))
							else:
								if not hasattr(self, "quiet"):
									self.warn(
										"{}{}{}{}".format("Warning: not both intron while parsing alt_cHGVS : ", tid,
														  ":", trannoEnt['c']))
						else:
							if not hasattr(self, "quiet"):
								self.warn(
									"{}{}{}{}".format("Warning: Unknown chgvs5 while parsing alt_cHGVS : ", tid, ": ",
													  chgvs_5))
						if renew_cHGVS is not None:
							if dupopt:
								trannoEnt['standard_cHGVS'] = renew_cHGVS
							else:
								trannoEnt['alt_cHGVS'] = renew_cHGVS

				if int(real_var.ref_cn) < int(real_var.alt_cn) and int(real_var.alt_cn) > 2:
					# non-duplication insertion

					ins_cn = int(real_var.alt_cn) - int(real_var.ref_cn)
					ins_cont = str(trRep) * int(ins_cn)

					if not dupopt:
						trannoEnt['standard_cHGVS'] = "{}{}{}{}{}{}".format(f, chgvs_5, trRep, '[', real_var.alt_cn,
																			']')

					m1 = re.search(r'^([\*\+]?)(\-?\d+)$', str(chgvs_3))
					m2 = re.search(r'^(\-?[^\-\+]+)(\+?d?)(\-?u?\d+)$', str(chgvs_3))
					if re.search(r'^\d+$', str(chgvs_3)):  # cds / ncRNA exon
						if (cdsOpt and int(chgvs_3) + 1 <= int(trdbEnt['csto']) - int(trdbEnt['csta'])) or (
							(not cdsOpt) and (chgvs_3 + 1 <= trdbEnt['len'])):
							trannoEnt['alt_cHGVS'] = "{}{}{}{}{}{}".format(f, chgvs_3, '_', int(chgvs_3) + 1, 'ins',
																		   ins_cont)
						else:
							outlatter = None
							if cdsOpt:
								outlatter = "*1"
							else:
								outlatter = "+1"
							trannoEnt['alt_cHGVS'] = "{}{}{}{}{}{}".format(f, chgvs_3, '_', outlatter, 'ins', ins_cont)
					elif m1:
						# utr or ncRNA ext
						sig = m1.group(1)
						ins_pos = m1.group(2)
						if str(ins_pos) == '-1':
							trannoEnt['alt_cHGVS'] = "{}{}{}".format(f, '-1_1ins', ins_cont)
						else:
							trannoEnt['alt_cHGVS'] = "{}{}{}{}{}{}{}".format(f, chgvs_3, '_', sig, int(ins_pos) + 1, 'ins',
																			 ins_cont)
					elif m2:
						# intron
						anchor = m2.group(1)
						sig = m2.group(2)
						ofst = m2.group(3)
						u3opt = 0
						if re.search(r'u', str(ofst)):
							u3opt = 1
							ofst = re.sub(r'u', '', str(ofst))
						if str(ofst) == '-1':
							trannoEnt['alt_cHGVS'] = "{}{}{}{}{}{}".format(f, chgvs_3, '_', anchor, 'ins', ins_cont)
						else:
							anc3 = int(ofst) + 1
							if u3opt:
								anc3 = re.sub(r'^-', '-u', str(anc3))
							trannoEnt['alt_cHGVS'] = "{}{}{}{}{}{}{}{}".format(f, chgvs_3, '_', anchor, sig, anc3,
																			   'ins', ins_cont)
					else:
						if not hasattr(self, 'quiet'):
							self.warn("Warning: Unknown chgvs5 while parsing alt_cHGVS : " + tid + ": " + chgvs_3)
			else:
				if int(cmpPos) == 0:  # 1 bp
					trannoEnt['c'] = str(f) + str(chgvs_5)
					if int(real_al) == 1:
						trannoEnt['c'] = "{}{}{}{}".format(trannoEnt['c'], real_r, '>', real_a)  # 1bp alt
					elif int(real_al) == 0:
						tmp = None
						if re.search(r'^[ACGTN]+$', str(real_r)):
							tmp = real_r
						else:
							tmp = ""
						trannoEnt['c'] = "{}{}{}".format(trannoEnt['c'], 'del', tmp)
					else:
						tmp = None
						if re.search(r'^[ACGTN]+$', str(real_r)):
							tmp = real_r
						else:
							tmp = ""
						trannoEnt['c'] = "{}{}{}{}{}".format(trannoEnt['c'], 'del', tmp, 'ins', real_a)
				elif int(cmpPos) > 0:  # multiple bp
					trannoEnt['c'] = "{}{}{}{}{}".format(f, chgvs_5, '_', chgvs_3, 'del')
					tmp = None
					if re.search(r'^[ACGTN]+$', str(real_r)) and len(real_r) < 50:
						tmp = real_r
					else:
						tmp = ''
					trannoEnt['c'] = "{}{}".format(trannoEnt['c'], tmp)
					if int(real_al) > 0:
						trannoEnt['c'] = "{}{}{}".format(trannoEnt['c'], 'ins', real_a)
				else:  # ins : reverse the positions
					trannoEnt['c'] = "{}{}{}{}{}{}".format(f, chgvs_3, '_', chgvs_5, 'ins', real_a)

			# 3. check if transcript-ablation
			if (re.search(r'^\-', str(trBegin)) or (str(trBegin) == '1')) and (re.search(r'^\+', str(trEnd)) or str(trEnd) == str(trdbEnt['len'])):
				trannoEnt['func'] = 'knockout'

				# need prot begin end?

				continue

			# * check if exon remapping introduced special case
			if re.search(r'^\+|\+d', str(chgvs_5)) and re.search(r'^\+|\+d', str(chgvs_3)):
				trannoEnt['func'] = 'unknown'
				continue

			if "preStart" in trannoEnt and "postEnd" in trannoEnt and str(trannoEnt['preStart']['exin']) != str(trannoEnt['postEnd']['exin']) and re.search(r'^IVS|^\.', str(trannoEnt['preStart']['exin'])) and (re.search(r'^IVS', str(trannoEnt['postEnd']['exin'])) or re.search(r'^3D', str(trannoEnt['postEnd']['r']))):
				# highest priority
				# due to the possibility of only delete one exon
				# without any intron region, which can be further
				# analysis the protein change, so we continue
				# without skip the following steps
				trannoEnt['alt_func'] = 'exon-loss'

			# 4. check if span different exon/intron/promoter/downstream
			if (re.search(r'^\+|\+d', str(chgvs_3)) or trannoEnt['ei_Begin'] != trannoEnt['ei_End']) and int(cmpPos) >= 0:  # cmpPos >= 0  not insertion at the edge
				trannoEnt['func'] = 'span'

				# need prot begin end?

				continue

			# 5. check if all promoter
			if str(trannoEnt['r_Begin']) == 'PROM' and trannoEnt['r_Begin'] == trannoEnt['r_End']:
				trannoEnt['func'] = 'promoter'
				continue

			# 6. check if all intron
			if re.search(r'IVS', str(trannoEnt['ei_Begin'])) and trannoEnt['ei_Begin'] == trannoEnt['ei_End']:
				if trannoEnt['genepartSO'] == 'abnormal-intron':
					trannoEnt['func'] = 'abnormal-intron'
				else:
					if hasattr(self, "genome_h") and re.search(r'^[AD]', str(trannoEnt['r_Begin'])) and trannoEnt['r_End'] == trannoEnt['r_Begin'] and len(trRef) == len(trAlt):
						# only in splice site substitution, to check if conanical
						cis_tag = trannoEnt['r_Begin'][:1]
						ext_len = 2 - len(trRef)
						if ext_len < 0:
							self.throw("Error length for trRef")
						Ladded, Radded = "", ""
						gchr, gstart, gend = annoEnt['var']['chr'], annoEnt['var']['pos'], annoEnt['var']['end']

						if not re.search(r'^chr', str(gchr)):
							gchr = "chr" + gchr

						if ext_len == 1:  # only can be 1
							if (trannoEnt['strd'] == '+' and (
								re.search(r'\+2$', trannoEnt['rnaBegin']) or re.search(r'\-1$',
																					   trannoEnt['rnaBegin']))) or (
									trannoEnt['strd'] == '-' and (
								re.search(r'\+1$', trannoEnt['rnaEnd']) or re.search(r'\-2$', trannoEnt['rnaEnd']))):
								extpos = gstart  # extend 1 bp left
								rgn_tmp = "{}{}{}{}{}".format(gchr, ":", extpos, "-", extpos)
								if rgn_tmp in getseq_cache:
									Ladded = getseq_cache[rgn_tmp]
								else:
									toAdd1 = None
									toAdd1 = self.genome_h.getseq(rgn_tmp)  # getseq函数找不到哪里实现的
									if toAdd1 is not None:
										Ladded = toAdd1.upper()
										getseq_cache[rgn_tmp] = Ladded
									else:
										self.warn(
											"{}{}".format("Warning: Can not get string from genome for ", rgn_tmp))
							else:
								extpos2 = gend + 1  # extend 1 bp right
								rgn_tmp2 = "{}{}{}{}{}".format(gchr, ":", extpos2, "-", extpos2)
								if rgn_tmp2 in getseq_cache:
									Radded = getseq_cache[rgn_tmp2]
								else:
									toAdd2 = None
									toAdd2 = self.genome_h.getseq(rgn_tmp2)  # getseq函数找不到哪里实现的
									if toAdd2 is not None:
										Radded = toAdd2.upper()
										getseq_cache[rgn_tmp2] = Radded
									else:
										self.warn(
											"{}{}".format("Warning: Can not get string from genome for ", rgn_tmp2))

							if trannoEnt['strd'] == '-':  # change to tr strand
								tmp = Ladded
								Ladded = self.rev_comp(Radded)
								Radded = self.rev_comp(tmp)

						toCheckRef, toCheckAlt = "{}{}{}".format(Ladded, trRef.upper(), Radded), "{}{}{}".format(Ladded,
																												 trAlt.upper(),
																												 Radded)

						if toCheckRef != BedAnno.canonicalSS["cis_tag"] and toCheckAlt == BedAnno.canonicalSS["cis_tag"]:
							trannoEnt['func'] = 'no-change'
							# here may need to keep the original format of hgvs.
							# $trannoEnt->{c}    = 'c.=';
							continue

					if re.search(r'^D', trannoEnt['r_Begin']) and re.search(r'^A', trannoEnt['r_End']):
						trannoEnt['func'] = 'splice'
					elif re.search(r'^D', trannoEnt['r_Begin']):
						trannoEnt['func'] = 'splice-5'
					elif re.search(r'^A', trannoEnt['r_End']):
						trannoEnt['func'] = 'splice-3'
					else:
						trannoEnt['func'] = 'intron'
				continue

			# 7. check if all exon or ex/intron edge insertion
			#    count all this insertion to be exon region insertion
			if (trannoEnt['ei_Begin'] == trannoEnt['ei_End'] and re.search(r'EX', trannoEnt['ei_Begin'])) or (cmpPos < 0 and trannoEnt['ei_Begin'] != trannoEnt['ei_End']):
				# any edge-insertion case will count on the exon change
				# any coding-utr edge insertion case will count on the
				# utr parts
				if not cdsOpt:
					if re.search(r'^\+', str(chgvs_5)):
						trannoEnt['func'] = 'unknown'
					elif re.search(r'^\-', str(chgvs_3)):
						trannoEnt['func'] = 'promoter'
					elif str(chgvs_5) == '1' and str(chgvs_3) == str(trdbEnt['len']):
						trannoEnt['func'] = 'knockout'
					else:
						trannoEnt['func'] = 'ncRNA'
					continue
				else:  # coding RNA
					m1 = re.search(r'^(\d+)\+1', str(chgvs_5))
					m2 = re.search(r'^\-(\d+)', str(chgvs_5))
					if re.search(r'\-u1$', str(chgvs_3)):
						trannoEnt['func'] = 'promoter'
						continue
					elif re.search(r'\+d1$', str(chgvs_5)):  # 3' downstream edge insertion
						trannoEnt['func'] = 'unknown'
					elif re.search(r'^\-', str(chgvs_3)):  # besides utr-5 - cds edge
						trannoEnt['func'] = 'utr-5'
						continue
					elif re.search(r'^\*', str(chgvs_5)):  # besides utr-3 - cds edge
						trannoEnt['func'] = 'utr-3'
						continue
					elif re.search(r'^\-', str(chgvs_5)) and re.search(r'^\*', str(chgvs_3)):
						utr5_len = None
						m = re.search(r'^\-(\d+)$', str(chgvs_5))
						if m:
							utr5_len = m.group(1)
						utr3_len = None
						m = re.search(r'^\-(\d+)$', str(chgvs_3))
						if m:
							utr3_len = m.group(1)
						if int(utr5_len) == int(trdbEnt['csta']) and int(utr3_len) == int(trdbEnt['len']) - int(trdbEnt['csto']):
							trannoEnt['func'] = 'knockout'
						else:
							trannoEnt['func'] = 'cds-loss'
						trannoEnt['p'] = 'p.0?'
						continue
					elif str(chgvs_3) == '1-1':
						# cds-begin - last 5'utr intron edge
						trannoEnt['func'] = 'utr-5'
						continue
					elif m1 and int(m1.group(1)) == int(trdbEnt['csto']) - int(trdbEnt['csta']):
						trannoEnt['func'] = 'utr-3'
						continue
					# here protein sequence var should be reparsed
					elif m2:
						u5_len = m2.group(1)
						if real_var.imp == 'rep':
							if int(u5_len) > int(real_rl) - int(real_al):
								trannoEnt['func'] = 'utr-5'
							else:
								trannoEnt['func'] = 'unknown'
						else:
							trannoEnt['func'] = 'init-loss'
							trannoEnt['P'] = 'p.0?'
						continue
					else:
						m = re.search(r'^\*(\d+)', str(chgvs_3))
						if m:
							# variant get across the 3' end of cds
							u3_len = m.group(1)
							if str(real_var.imp) == 'rep' and int(u3_len) > (int(real_rl) - int(real_al)):
								trannoEnt['func'] = 'utr-3'
								# hgvs 2.1511 don't refer to this kind of
								# coding synon, just keep it.
								trannoEnt['p'] = 'p.(=)'
								continue
						m = re.search(r'^(\d+)\+1$', str(chgvs_5))
						if int(cmpPos) < 0 and m:
							chgvs_5 = int(m.group(1)) + 1
						m = re.search(r'^(\d+)\-1$', str(chgvs_3))
						if int(cmpPos) < 0 and m:
							chgvs_3 = int(m.group(1)) - 1
						m = re.search(r'^(\d+)\+1$', str(trBegin))
						if int(cmpPos) < 0 and m:
							trBegin = int(m.group(1)) + 1
						m = re.search(r'^(\d+)\-1$', str(trEnd))
						if int(cmpPos) < 0 and m:
							trEnd = int(m.group(1)) - 1

						# 5' end should be in coding region.
						# leave this as debug information
						if not re.search(r'^\d+$', str(chgvs_5)) or int(chgvs_5) == 0:
							self.throw("{}{}".format("Error: unavailable 5' cHGVS. ", chgvs_5))
						rcInfo5, rcInfo3 = BedAnno._getPairedCodon(trdbEnt, trBegin, trEnd)
						if str(real_a) == "":
							# for deletion case try to predict the pr change
							if 5 < len(rcInfo5):
								rcInfo5.pop(4)
							if 5 < len(rcInfo3):
								rcInfo3.pop(4)

							if rcInfo5[0] > 0 and rcInfo3[0] > 0 and rcInfo5[4] == -1 and rcInfo3[4] == -1:
								hit_whole_deleted_fs = 0
								for fsld in trdbEnt['nfs'].keys():
									if int(fsld) == int(trBegin) - 1 and int(fsld) + int(trdbEnt['nfs'][fsld]) == int(trEnd):
										hit_whole_deleted_fs = 1
										break
								if hit_whole_deleted_fs:
									trannoEnt['protBegin'] = rcInfo5[0]
									trannoEnt['protEnd'] = rcInfo3[0]
									trannoEnt['func'] = 'no-change'
									continue

						if rcInfo5[0] > 0 and rcInfo5[4] == -1 and rcInfo3[0] > 0 and rcInfo3[4] == -1:
							trannoEnt['protBegin'] = rcInfo5[0]
							trannoEnt['protEnd'] = rcInfo3[0]
							trannoEnt['func'] = 'abnormal-fs-site'
							continue

						real_rl_hidden_ofst = 0
						if "nfs" in trdbEnt:
							for fsld in trdbEnt['nfs'].keys():
								fsEd = int(fsld) + int(trdbEnt['nfs'][fsld])
								if (trdbEnt['nfs'][fsld] > 0 and fsld >= trBegin and fsEd <= trEnd) or (
										trdbEnt['nfs'][fsld] < 0 and (
									(fsEd >= trBegin and fsld < trEnd) or (fsld >= trEnd and fsEd < trBegin))):
									# contain fs site or in duplicated fs site
									real_rl_hidden_ofst =int(real_rl_hidden_ofst)+ int(trdbEnt['nfs'][fsld])
								else:
									if trdbEnt['nfs'][fsld] > 0 and fsld < trBegin and trBegin <= fsEd:
										# begin in deleted fs site
										real_rl_hidden_ofst += fsEd - trBegin
									if trdbEnt['nfs'][fsld] > 0 and fsld < trEnd and trEnd <= fsEd:
										# end in deleted fs site
										real_rl_hidden_ofst += trEnd - fsld
									if trdbEnt['nfs'][fsld] < 0 and (
										(fsld >= trBegin and fsEd < trBegin) or (fsld >= trEnd and fsEd < trEnd)):
										# start or end in duplicated fs site
										real_rl_hidden_ofst += trdbEnt['nfs'][fsld]

						translate_opts = dict()
						m = re.search(r'^NM_MT-', tid)
						if m:
							translate_opts['mito'] = 1
						if "A" in trdbEnt:
							translate_opts['polyA'] = 1
						if "X" in trdbEnt or "U" in trdbEnt:
							translate_opts['nostop'] = 1
						altcodon_opts = copy.deepcopy(translate_opts)
						if "nostop" in altcodon_opts:
							del altcodon_opts['nostop']
						if "pseq" not in trdbEnt:
							pseq, frame_next = BedAnno.translate(trdbEnt['cseq'], translate_opts)
							trdbEnt['pseq'] = pseq

						prBegin, prEnd, prRef, prAlt = None, None, None, None

						prBegin = rcInfo5[0]

						ready_to_add_5 = rcInfo5[1][:rcInfo5[4]]

						diff_ra = real_rl - real_al - real_rl_hidden_ofst

						# probably frame shift flag
						frameshift_flag = 1 if diff_ra % 3 > 0 else 0

						# end in cds or not.
						end_in_cds_flag = 0 if re.search(r'^\*', str(chgvs_3)) or int(chgvs_3) > int(int(trdbEnt['csto']) - 3) else 1

						ready_to_add_3 = None
						if not end_in_cds_flag or frameshift_flag:
							ready_to_add_3 = BedAnno._cdsubstr(trdbEnt, trEnd)
							ready_to_add_3_len =len(ready_to_add_3)
							# this variants's effect will be end at the terminal
							prEnd = int(trdbEnt['plen']) + 1
							if "A" in trdbEnt:
								ready_to_add_3 = "{}{}".format(ready_to_add_3, "A" * (3 - (len(trdbEnt['cseq']) % 3)))
						else:  # chgvs_3 in cds region with no frameshift
							ready_to_add_3 = rcInfo3[1][(rcInfo3[4] + 1):3]
							prEnd = rcInfo3[0]

						if int(trdbEnt['plen']) + 1 != len(trdbEnt['pseq']):
							if not hasattr(self, 'quiet'):
								self.warn("[Critical Warning] : annotation for " + tid + " may not be correct, \n" + \
										  "                     due to a mis-translated base in its \n" + \
										  "                     cds region, usually a frameshift on \n" + \
										  "                     TGA to GAX, please ignore any annotation \n" + \
										  "                     on latter part behind that frameshift.\n" + \
										  "                     protein annotation for " + str(
									tid) + " will be skipped\n" + \
										  "                     and the function of it will be annotation-fail\n" + \
										  "                     until the problem is fixed.")
							trannoEnt['func'] = 'annotation-fail'
							continue

						# Make protBegin and protEnd be coordinate to the
						# coordinates of nucl variants, extend the
						# frameshift's effect to the end of protein
						trannoEnt['protBegin'] = prBegin
						trannoEnt['protEnd'] = prEnd
						# make the order of func assignment as follow:
						# * total
						# 1. init-loss
						# 2. frameshift
						# 3. coding-synon
						#
						# * 1bp aa change
						# 0. abnormal-inseq-stop
						# 1. stop-retained
						# 2. stop-loss
						# 3. nonsense
						# 4. missense
						#
						# * other variants
						# 0. abnormal-inseq-stop
						# 1. stop-retained
						# 2. frameshift
						# 3. stop-loss
						# 4. stop-gain
						# 5. cds-ins
						# 6. cds-del
						# 7. cds-indel

						# extend the altered sequence to codon coordinates
						codon_alt = "{}{}{}".format(ready_to_add_5, real_a, ready_to_add_3)
						codon_alt_len=len(codon_alt)
						# here only issues variant start/stop on exon region
						# and 5' start position on the coding region
						# init codon flag
						hit_init_flag = 1 if rcInfo5[0] == 1 else 0

						# stop codon flag
						hit_stop_flag = 1 if trEnd > trdbEnt['csta'] + 3 * trdbEnt['plen'] else 0

						# check if an alterstart
						start_codons = {'ATG': -1}
						if "altstart" in trdbEnt:
							for k in trdbEnt['altstart'].keys():
								start_codons[k] = 1
						init_synon = 0  # indicate start codon search result
						if hit_init_flag:
							# check if altered sequence have new start codon
							for startCodon in sorted(start_codons.keys()):
								if codon_alt[:3] == startCodon:
									init_synon = 1
									break

							if init_synon == 0:
								if re.search(r'N', trAlt):
									all_hit = 1
									for base in ['A', 'C', 'G', 'T']:
										subalt_codon = codon_alt
										subalt_codon = re.sub(r'N', base, subalt_codon)
										if subalt_codon not in start_codons:
											all_hit = 0
											break

									if 3 == len(codon_alt):
										trannoEnt['cc'] = "{}{}{}".format(rcInfo5[1], "=>", codon_alt)
									if all_hit == 1:
										trannoEnt['func'] = "altstart"
										trannoEnt['prAlt'] = "M"
										init_synon = 1
									else:
										trannoEnt['func'] = 'unknown-no-call'
								else:
									trannoEnt['func'] = 'init-loss'
									trannoEnt['p'] = 'p.0?'
									if str(prBegin) == '1' and str(prEnd) == '1':
										zero = None
										trannoEnt['prRef'] = 'M'
										altCodon = codon_alt[:3]
										trannoEnt['prAlt'], zero = BedAnno.translate(altCodon, altcodon_opts)
									continue
							else:
								# can it continue to be assigned
								# other function code and pHGVS
								# or just finished as altstart?
								#
								# I choose to continue.
								trannoEnt['func'] = 'altstart'

						# give out the first of single codon to be
						codon_to_be = codon_alt[:3]

						if prBegin - 1 == trdbEnt['plen'] and prEnd == prBegin:
							prRef = '*'
						elif prBegin - 1 < trdbEnt['plen']:
							prRef = trdbEnt['pseq'][(prBegin - 1):(prEnd)]
						else:
							self.throw("".join(map(str, ["Error: Cannot get prRef from " + tid + " [" + trdbEnt[
								'plen'] + prBegin + prEnd + "]\n" + \
														 "      Var: [" + annoEnt['var']['chr'] + annoEnt['var'][
															 'pos'] + annoEnt['var']['end'] + \
														 annoEnt['var']['ref'] + annoEnt['var']['alt']])))

						# we can not allow inseq stop codon in altered sequence
						# due to frameshift and ambiguity.

						next_alt_frame = None
						prAlt, next_alt_frame = BedAnno.translate(codon_alt, altcodon_opts)

						if hit_init_flag:
							prRef = re.sub(r'^[A-Z]', 'M', prRef)  # change init pep to M      ,change alt-init to M
							if init_synon:
								prAlt = re.sub(r'^[A-Z]', 'M', prAlt)

						trannoEnt['prRef'] = prRef
						if 'prAlt' not in trannoEnt:
							trannoEnt['prAlt'] = prAlt

						prAlt = trannoEnt['prAlt']  # uniform prAlt

						# to indicate whether the alternate sequence
						# encode a stop codon or non-frameshift, or otherwise
						# with a non stopped frameshift.
						non_stop_flag = 1 if (not re.search(r'\*$', prAlt)) and (
						prRef == "*" or next_alt_frame or frameshift_flag) else 0
						# To avoid large range perlre, before parse protein
						# var we should first deal with frameshift issues
						# need to deal with non_stop_flag and frameshift and
						# imp ref case.
						#
						# assign cc and polar to the same length,
						# single aa substitution, no matter which position.
						if diff_ra == 0 and rcInfo5[0] == rcInfo3[0]:
							# using trim due to terminal codon
							aa_to_be, polar_to_be = BedAnno._getAAandPolar(codon_to_be, altcodon_opts)

							trannoEnt['cc'] = "{}{}{}".format(rcInfo5[1], '=>', codon_to_be)
							if prRef != prAlt:
								trannoEnt['polar'] = "{}{}{}".format(rcInfo5[3], '=>', polar_to_be)

						prSimpleSame = 0
						pp = 0
						while pp < len(prRef) and pp < len(prAlt):
							if prRef[pp:(pp + 1)] == prAlt[pp:(pp + 1)]:
								prSimpleSame += 1
							else:
								break
							pp += 1
						no_parsed_pP = prBegin - 1 + prSimpleSame
						no_parsed_prStart = trdbEnt['pseq'][no_parsed_pP:(no_parsed_pP + 1)]
						no_parsed_prAlt_Start = prAlt[prSimpleSame:(prSimpleSame + 1)]

						# non stop frameshift with extra non-coded base
						if non_stop_flag:
							if int(no_parsed_pP) >= int(trdbEnt['plen']):
								trannoEnt['func'] = 'stop-loss'
							else:
								trannoEnt['func'] = 'frameshift'
							trannoEnt['p'] = "{}{}{}{}{}".format('p.', no_parsed_prStart, int(no_parsed_pP) + 1,
																 no_parsed_prAlt_Start, 'fs*?')
							continue

						# identical to reference (can be from any kind of vars)
						if prRef == prAlt:
							if init_synon:  # altstart
								# do nothing
								pass
							elif hit_stop_flag:
								# need this special assertion?
								trannoEnt['func'] = 'stop-retained'
							else:
								trannoEnt['func'] = 'coding-synon'
							if 1 == len(prRef):
								trannoEnt['p'] = "{}{}{}{}".format('p.(', prRef, no_parsed_pP, '=)')
							else:
								trannoEnt['p'] = "{}{}{}{}{}".format('p.(', no_parsed_pP, '_',
																	 (no_parsed_pP + len(prRef) - 1), '=)')
							continue

						# frameshift
						if not end_in_cds_flag or frameshift_flag:
							trannoEnt['func'] = 'frameshift'
							trannoEnt['p'] = "{}{}{}{}".format('p.', no_parsed_prStart, no_parsed_pP + 1,
															   no_parsed_prAlt_Start)
							ext_length = len(prAlt) - prSimpleSame
							if ext_length > 1 or not re.search(r'\*$', prAlt):
								trannoEnt['p'] += 'fs*'
								if re.search(r'\*$', prAlt):  # ext length estimated
									trannoEnt['p'] += str(ext_length)
								else:  # don't meet a stop codon
									trannoEnt['p'] += '?'
							trannoEnt['protBegin'] += int(prSimpleSame)
							continue

						ins_stop_tag = 0
						if not re.search(r'\*', prRef) and prAlt != '*' and re.search(r'\*$', prAlt):
							# stop gain
							prEnd = trdbEnt['plen']  # extent to the end of prot
							prRef = trdbEnt['pseq'][(int(prBegin) - 1):int(prEnd)]
							prAlt = re.sub(r'\*$', '', prAlt)
							ins_stop_tag = 1
						if prAlt == '?':
							trannoEnt['func'] = 'unknown-no-call'
							trannoEnt['p'] = 'p.'
							if prBegin == prEnd:
								trannoEnt['p'] = "{}{}{}{}".format(trannoEnt['p'], prRef, prBegin, '?')
							elif prBegin < prEnd:
								trannoEnt['p'] = "{}{}{}{}{}{}{}".format(trannoEnt['p'], prRef[:1], prBegin, '_',
																		 prRef[(-1):], prEnd, 'delins?')
							else:
								trannoEnt['p'] = "{}{}{}{}{}{}{}".format(trannoEnt['p'],
																		 trdbEnt['pseq'][(prEnd - 1):prEnd], prEnd, '_',
																		 trdbEnt['pseq'][(prBegin - 1):prBegin],
																		 prBegin, 'ins?')
							continue

						# parse the protein variants
						# to recognize the repeat and adjust to correct position
						prVar = self.prWalker(trdbEnt['pseq'], prBegin, prEnd, prRef, prAlt)

						# 0-based
						p_P, p_r, p_a, prl, pal = prVar.getUnifiedVar('+')
						prStart = trdbEnt['pseq'][p_P:(p_P + 1)]
						prStop = None
						if p_P > 0 or prl > 0:
							prStop = trdbEnt['pseq'][(p_P + prl - 1):(p_P + prl)]

						# single substitution
						if prVar.imp == 'snv':
							if p_r == 'X':
								trannoEnt['func'] = 'abnormal-inseq-stop'
							elif p_r == '.':  # N in refseq transcript
								trannoEnt['func'] = 'unknown'
							elif p_a == '?':  # substitution with N
								trannoEnt['func'] = 'unknown-no-call'
							elif p_r == '*':
								trannoEnt['func'] = 'stop-loss'
							elif p_a == '*':
								trannoEnt['func'] = 'nonsense'
							elif init_synon:
								trannoEnt['func'] = 'altstart'
							else:
								trannoEnt['func'] = 'missense'
							trannoEnt['p'] = "{}{}{}{}".format('p.', p_r, (p_P + 1), p_a)
							continue

						# repeat
						if prVar.imp == 'rep':
							if prVar.ref_cn < prVar.alt_cn:
								trannoEnt['func'] = 'cds-ins'
							else:
								trannoEnt['func'] = 'cds-del'
							if prVar.replen == 1:
								trannoEnt['p'] = "{}{}{}".format('p.', prVar.rep, p_P + 1)
							else:
								trannoEnt['p'] = "{}{}{}{}{}{}".format('p.', prStart, p_P + 1, '_', prVar.rep[(-1):],
																	   p_P + prVar.replen)
							if prVar.ref_cn == 1 and prVar.alt_cn == 2:
								trannoEnt['p'] = "{}{}".format(trannoEnt['p'], 'dup')
							else:
								trannoEnt['p'] = "{}{}{}{}{}{}".format(trannoEnt['p'], '[', prVar.ref_cn, '>',
																	   prVar.alt_cn, ']')

							# add a new key alt_pHGVS for querying
							# previous database, do exactly like alt_cHGVS
							# but will be much more simple.

							phgvs_5 = p_P + 1
							phgvs_3 = p_P + (prVar.ref_cn * prVar.replen)
							repPrSta = prVar.rep[:1]
							repPrEnd = prVar.rep[-1:]
							if prVar.ref_cn > prVar.alt_cn:
								# deletion or duplication in previous definition
								if prVar.ref_cn - prVar.alt_cn == 1 and prVar.replen == 1:
									trannoEnt['alt_pHGVS'] = "{}{}{}{}".format('p.', prVar.rep, phgvs_3, 'del')
								else:
									trannoEnt['alt_pHGVS'] = "{}{}{}{}{}{}{}".format('p.', repPrSta,
																					 phgvs_5 + prVar.alt_cn * prVar.replen,
																					 '_', repPrEnd, phgvs_3, 'del')
							else:  # insertion
								ins_cn = prVar.alt_cn - prVar.ref_cn
								ins_cont = prVar.rep * ins_cn  # 注意prVar.rep是数字还是字符串，*对他们的结果是不同的
								prPost = trdbEnt['pseq'][phgvs_3:phgvs_3 + 1]
								trannoEnt['alt_pHGVS'] = "{}{}{}{}{}{}{}{}".format('p.', repPrEnd, phgvs_3, '_', prPost,
																				   phgvs_3 + 1, 'ins', ins_cont)
							continue

						# when inserted bases contains stop codon but not framshift
						if re.search(r'^\*$', p_a):  # prAlt is only a stop codon
							trannoEnt['p'] = "{}{}{}{}".format('p.', prStart, p_P + 1, '*')
							trannoEnt['func'] = 'stop-gain'
							continue

						# insertion
						if prVar.sm == 0:

							# insert into the 5'edge of start codon
							# this may caused by a insertion or delins
							# around init-codon to recreat a new init-codon
							# in the altered sequence
							if p_a == '?':  # substitution with N
								trannoEnt['func'] = 'unknown-no-call'
							elif p_P == 0:

								# altstart don't fix this
								# use pHGVS to indicate
								# altstart here
								#
								# hgvs 2.1511 don't refer to this kind of
								# coding synon, so just keep it.
								trannoEnt['p'] = "p.(=)"

								# can it be utr-5?
								trannoEnt['func'] = 'utr-5'
								continue
							# insert between the last aa and stop codon
							# this may caused by a insertion or delins
							# aroud the stop codon
							elif hit_stop_flag:
								trannoEnt['func'] = 'stop-retained'
							else:
								trannoEnt['func'] = 'stop-gain' if ins_stop_tag else 'cds-ins'

							trannoEnt['p'] = "{}{}{}{}{}{}{}{}".format('p.', prStop, p_P, '_', prStart, p_P + 1, 'ins',
																	   p_a)
							continue

						if prVar.sm >= 1:
							if p_a == '?':  # substitution with N
								trannoEnt['func'] = 'unknown-no-call'
							elif re.search(r'\*', p_r):
								trannoEnt['func'] = 'stop-loss'
							elif hit_stop_flag:
								trannoEnt['func'] = 'stop-retained'
							elif ins_stop_tag:
								trannoEnt['func'] = 'stop-gain'
							elif pal == 0:
								trannoEnt['func'] = 'cds-del'
							elif pal == prl:
								trannoEnt['func'] = 'missense'
							else:
								trannoEnt['func'] = 'cds-indel'
							trannoEnt['p'] = "{}{}{}".format('p.', prStart, p_P + 1)
							if prVar.sm == 1:
								trannoEnt['p'] = "{}{}".format(trannoEnt['p'], 'del')
							else:
								trannoEnt['p'] = "{}{}{}{}{}".format(trannoEnt['p'], '_', prStop, p_P + prl, 'del')
							if p_a != "":
								trannoEnt['p'] = "{}{}{}".format(trannoEnt['p'], 'ins', p_a)
							continue
							# any other cases?

		return annoEnt

	'''
	=head2 finaliseAnno

		About   : finalise the BedAnno::Anno entry by check all tag values,
				  and uniform them for AE output usage, query transcript
				  oringinated additional resources and add them into the data
				  frame.
		Usage   : $beda->finaliseAnno($annEnt);
		Args    : BedAnno entry and a BedAnno::Anno entry
		Returns : A finalised BedAnno::Anno entry

	=cut
	'''
	def finaliseAnno(self, annoEnt):
		#if "trInfo" in annoEnt:
		if hasattr(annoEnt,"trInfo"):

			# use this extended local db to check if exon variant locates
			# in the 3bp range of exon-intron edge.
			av = annoEnt.var
			tmp = None
			if int(av.pos) - 3 > 0:
				tmp = int(av.pos) - 3
			else:
				tmp = 0
			ext_region = "{}{}{}{}{}".format(av.chr, ':', tmp, '-', int(av.pos) + int(av.reflen) + 3)
			ext_localdb = self.load_anno({"region": ext_region})
			ext_localdb = ext_localdb[av.chr]  # no need to check exists

			for tid in sorted(annoEnt.trInfo.keys()):

				trAnnoEnt = annoEnt.trInfo[tid]

				if re.search(r'^N[MR]_MT-', tid):
					trAnnoEnt['geneSym'] = re.sub(r'^MT-', '', trAnnoEnt['geneSym'])

				qtid = tid
				qtid = re.sub(r'\-\d+$', '', qtid)
				trdbEnt = self.trInfodb[qtid]
				trAnnoEnt['prot'] = "" if ("prot" not in trdbEnt or trdbEnt['prot'] == '.') else trdbEnt['prot']

				strd = 1 if trAnnoEnt['strd'] == '+' else 0

				# complete informations depend on already assigned values
				genepartSO = None
				if str(trAnnoEnt['rnaBegin']) == '?' or str(trAnnoEnt['rnaEnd']) == '?' or str(trAnnoEnt['genepartSO']) == "annotation-fail":
					# do nothing?
					genepartSO = trAnnoEnt['genepartSO']
					trAnnoEnt['r'] = '?'
					trAnnoEnt['exin'] = '?'
					trAnnoEnt['componentIndex'] = ''
					m = re.search(r'EX(\d+)', trAnnoEnt['ei_Begin'])
					if m:
						trAnnoEnt['exonIndex'] = m.group(1)
					m = re.search(r'IVS(\d+)', trAnnoEnt['ei_Begin'])
					if m:
						trAnnoEnt['intronIndex'] = m.group(1)
				elif trAnnoEnt['r_Begin'] == trAnnoEnt['r_End']:
					if re.search(r'^\d+$', str(trAnnoEnt['genepartSO'])):
						genepartSO = "SO:%07d" % int(trAnnoEnt['genepartSO'])
					else:
						genepartSO = trAnnoEnt['genepartSO']
					if "r" not in trAnnoEnt:
						trAnnoEnt['r'] = trAnnoEnt['r_Begin']
					if "exin" not in trAnnoEnt:
						trAnnoEnt['exin'] = trAnnoEnt['ei_Begin']
					m1 = re.search(r'^EX(\d+)', trAnnoEnt['exin'])
					m2 = re.search(r'^IVS(\d+)', trAnnoEnt['exin'])
					if m1:
						trAnnoEnt['componentIndex'] = m1.group(1)
						trAnnoEnt['exonIndex'] = m1.group(1)
					elif m2:
						trAnnoEnt['componentIndex'] = m2.group(1)
						trAnnoEnt['intronIndex'] = m2.group(1)
					else:  # PROM
						trAnnoEnt['componentIndex'] = 0
				elif re.search(r'^IVS', trAnnoEnt['ei_Begin']) and trAnnoEnt['ei_Begin'] == trAnnoEnt[
					'ei_End']:  # span splice or interior intron
					if re.search(r'^D', trAnnoEnt['r_Begin']) and re.search(r'^A', trAnnoEnt['r_End']):
						genepartSO = 'span'
					elif re.search(r'^D', trAnnoEnt['r_Begin']):
						genepartSO = "SO:0000163"  # 5'
					elif re.search(r'^A', trAnnoEnt['r_End']):
						genepartSO = "SO:0000164"  # 3'
					elif re.search(r'^A', trAnnoEnt['r_Begin']):  # insertion on intron edge
						genepartSO = BedAnno._intronSO(trAnnoEnt['r_End'])
					elif re.search(r'^D', trAnnoEnt['r_End']):
						genepartSO = BedAnno._intronSO(trAnnoEnt['r_Begin'])
					else:
						raise "Unknown Error in span case intron!"
					trAnnoEnt['r'] = "{}{}{}".format(trAnnoEnt['r_Begin'], '-', trAnnoEnt['r_End'])
					trAnnoEnt['exin'] = trAnnoEnt['ei_Begin']
					m = re.search(r'(\d+)', trAnnoEnt['ei_Begin'])
					if m:
						trAnnoEnt['componentIndex'] = m.group(1)
					trAnnoEnt['intronIndex'] = trAnnoEnt['componentIndex']
				else:
					genepartSO = 'span'
					trAnnoEnt['componentIndex'] = ''
					m = re.search(r'(\d+)', trAnnoEnt['ei_Begin'])
					if m:  # index the begin part
						trAnnoEnt['componentIndex'] = m.group(1)
					trAnnoEnt['r'] = "{}{}{}".format(trAnnoEnt['r_Begin'], '-', trAnnoEnt['r_End'])
					if trAnnoEnt['ei_Begin'] == trAnnoEnt['ei_End']:
						trAnnoEnt['exin'] = trAnnoEnt['ei_Begin']
						if trAnnoEnt['componentIndex'] != '':
							trAnnoEnt['exonIndex'] = trAnnoEnt['componentIndex']
					else:
						trAnnoEnt['exin'] = "{}{}{}".format(trAnnoEnt['ei_Begin'], '-', trAnnoEnt['ei_End'])
						m = re.search(r'EX(\d+)', trAnnoEnt['ei_Begin'])
						if m:
							trAnnoEnt['exonIndex'] = m.group(1)
						m = re.search(r'IVS(\d+)', trAnnoEnt['ei_Begin'])
						if m:
							trAnnoEnt['intronIndex'] = m.group(1)
						if trAnnoEnt['r_Begin'] == 'PROM' and re.search(r'^5U', trAnnoEnt['r_End']):
							trAnnoEnt['componentIndex'] = 0
						else:  # non-equal exon intron number or UTR/CDS span
							if trAnnoEnt['trRef'] == '':  # edge insertion case
								m = re.search(r'EX(\d+)', trAnnoEnt['exin'])
								if m:
									trAnnoEnt['exonIndex'] = m.group(1)
								m = re.search(r'IVS(\d+)', trAnnoEnt['exin'])
								if m:
									trAnnoEnt['intronIndex'] = m.group(1)

				if genepartSO not in BedAnno.SO2Name:
					raise "Error: unknown genepartSO [{}].".format(genepartSO)

				if "exin" not in trAnnoEnt:
					trAnnoEnt['exin'] = '.'
				if "r" not in trAnnoEnt:
					trAnnoEnt['r'] = '.'
				trAnnoEnt['genepart'] = BedAnno.SO2Name[genepartSO]
				trAnnoEnt['genepartSO'] = genepartSO if re.search(r'^SO:', genepartSO) else ""
				if "componentIndex" not in trAnnoEnt:
					trAnnoEnt['componentIndex'] = ''
				if "exonIndex" not in trAnnoEnt:
					trAnnoEnt['exonIndex'] = '.'
				if "intronIndex" not in trAnnoEnt:
					trAnnoEnt['intronIndex'] = '.'

				# change annotation-fail 's other value to empty?

				# uniform protBegin end
				if "protBegin" not in trAnnoEnt or trAnnoEnt['protBegin'] == '0':
					trAnnoEnt['protBegin'] = ""
				else:
					trAnnoEnt['protBegin'] = trAnnoEnt['protBegin']
				if "protEnd" not in trAnnoEnt or trAnnoEnt['protEnd'] == '0':
					trAnnoEnt['protEnd'] = ""
				else:
					trAnnoEnt['protEnd'] = trAnnoEnt['protEnd']

				if "func" not in trAnnoEnt:
					trAnnoEnt['func'] = 'unknown'

				# uniform function group
				if trAnnoEnt['func'] not in BedAnno.func2SO:
					raise "Error: unknown func code [{}].".format(trAnnoEnt['func'])
				trAnnoEnt['funcSO'] = BedAnno.func2SO[trAnnoEnt['func']]
				trAnnoEnt['funcSOname'] = BedAnno.SO2Name[trAnnoEnt['funcSO']]
				trAnnoEnt['funcSO'] = trAnnoEnt['funcSO'] if re.search(r'^SO:', trAnnoEnt['funcSO']) else ""
				if "c" in trAnnoEnt and (
					re.search(r'\d+\+5[ACGT]>[ACGT\?]$', trAnnoEnt['c']) or re.search(r'\d+\+5del[ACGT]?$',
																					  trAnnoEnt['c'])):
					# only single position
					trAnnoEnt['alt_func'] = 'splice-5-5th'
				elif not re.search(r'splice', trAnnoEnt['func']) and trAnnoEnt[
					'genepart'] != "annotation-fail" and "alt_func" not in trAnnoEnt:
					if trAnnoEnt['func'] == 'span' and (
								re.search(r'^IVS', trAnnoEnt['ei_Begin']) or re.search(r'^IVS', trAnnoEnt['ei_End']) or (
								trAnnoEnt['r_Begin'] == 'PROM' and not re.search(r'^EX1E?$', trAnnoEnt['ei_End'])) or (
							trAnnoEnt['r_End'] == '3D' and not re.search(r'^EX(\d+)E$', trAnnoEnt['ei_Begin']))):
						trAnnoEnt['alt_func'] = 'splice-region'  # or just ommit this?
					else:
						m = re.search(r'\d+[\+\-](\d+)', str(trAnnoEnt['rnaBegin']))
						if m:
							if int(m.group(1)) <= 8:
								trAnnoEnt['alt_func'] = 'splice-region'
							elif int(m.group(1)) <= 10:
								trAnnoEnt['alt_func'] = 'splice-ext'
						m = re.search(r'\d+[\+\-](\d+)', str(trAnnoEnt['rnaEnd']))
						if m:
							if int(m.group(1)) <= 8:
								trAnnoEnt['alt_func'] = 'splice-region'
							elif int(m.group(1)) <= 10:
								if "alt_func" not in trAnnoEnt:
									trAnnoEnt['alt_func'] = 'splice-ext'

						if "alt_func" not in trAnnoEnt:  # on exon region
							for k in range(len(ext_localdb)):
								if "detail" not in ext_localdb[k]:
									ext_localdb[k] = self.assign_detail(ext_localdb[k])
								if tid not in ext_localdb[k]['detail']:
									continue
								rtd = ext_localdb[k]['detail'][tid]
								if re.search(r'^[AD]', rtd['blka']):  # hit splice site
									# calculate distance from exon-intron edge
									if BedAnno._checkDistEdge(trAnnoEnt, rtd):
										trAnnoEnt['alt_func'] = 'splice-region'
										break

				if "alt_func" in trAnnoEnt:
					trAnnoEnt['alt_funcSO'] = BedAnno.func2SO[trAnnoEnt['alt_func']]
					trAnnoEnt['alt_funcSOname'] = BedAnno.SO2Name[trAnnoEnt['alt_funcSO']]
					trAnnoEnt['alt_funcSO'] = trAnnoEnt['alt_funcSO'] if re.search(r'^SO:', trAnnoEnt['alt_funcSO']) else ""

				# add additional resource
				if "prot" in trAnnoEnt and str(trAnnoEnt['prot']) != "":
					trAnnoEnt['protBegin'] = "" if str(trAnnoEnt['protBegin']) == '0' else trAnnoEnt['protBegin']
					trAnnoEnt['protEnd'] = "" if str(trAnnoEnt['protEnd']) == '0' else trAnnoEnt['protEnd']
					if str(trAnnoEnt['protBegin']) != "" and str(trAnnoEnt['protEnd']) != "":
						if hasattr(self, "pfam"):
							pb, pe = trAnnoEnt['protBegin'], trAnnoEnt['protEnd']
							if int(pb) > int(pe):
								tmp = pb
								pb = pe
								pe = tmp
							trAnnoEnt['pfamId'], trAnnoEnt['pfamName'] = self.pfam_h.getPfam(trAnnoEnt['prot'], pb,pe)  # getPfam函数的实现不知道在哪里
					if hasattr(self, "prediction"):
						m = re.search(r'^p\.[A-Z](\d+)([A-Z])$', str(trAnnoEnt['p']))
						if "p" in trAnnoEnt and m:
							rpred = self.prediction_h.getPredScore(trAnnoEnt['prot'], m.group(1),m.group(2))  # getPredScore函数的实现不知道在哪里
							if "sift" in rpred:
								trAnnoEnt['siftPred'] = rpred['sift'][0]
								trAnnoEnt['siftScore'] = rpred['sift'][1]
							if "polyphen2_humdiv" in rpred:
								trAnnoEnt['pp2divPred'] = rpred['polyphen2_humdiv'][0]
								trAnnoEnt['pp2divScore'] = rpred['polyphen2_humdiv'][1]
							if "polyphen2_humvar" in rpred:
								trAnnoEnt['pp2varPred'] = rpred['polyphen2_humvar'][0]
								trAnnoEnt['pp2varScore'] = rpred['polyphen2_humvar'][1]
						elif "prRef" in trAnnoEnt and "prAlt" in trAnnoEnt:
							prReflen = len(trAnnoEnt['prRef'])
							prAltlen = len(trAnnoEnt['prAlt'])
							if prReflen == 1 and prAltlen == 1 and trAnnoEnt['protBegin'] == '1' and trAnnoEnt[
								'protEnd'] == '1':
								init_pred = self.prediction_h.getPredScore(trAnnoEnt['prot'], 1, trAnnoEnt['prAlt'])
								if "sift" in init_pred:
									trAnnoEnt['siftPred'] = rpred['sift'][0]
									trAnnoEnt['siftScore'] = rpred['sift'][1]
								if "polyphen2_humdiv" in rpred:
									trAnnoEnt['pp2divPred'] = rpred['polyphen2_humdiv'][0]
									trAnnoEnt['pp2divScore'] = rpred['polyphen2_humdiv'][1]
								if "polyphen2_humvar" in rpred:
									trAnnoEnt['pp2varPred'] = rpred['polyphen2_humvar'][0]
									trAnnoEnt['pp2varScore'] = rpred['polyphen2_humvar'][1]

						if hasattr(self, "condel_h") and self.condel_h:
							rcondel_info = self.condel_h.pred(trAnnoEnt)  # pred函数的实现不知道在哪里
							if "condelPred" in rcondel_info and rcondel_info['condelPred'] != "not_computable_was":
								trAnnoEnt['condelPred'] = rcondel_info['condelPred']
								trAnnoEnt['condelScore'] = rcondel_info['condelScore']

				if "p" in trAnnoEnt and trAnnoEnt['p'] and re.search(r'\]$', trAnnoEnt['p']):
					trAnnoEnt['standard_pHGVS'] = BedAnno.gen_standard_gphgvs(trAnnoEnt['p'])

				if "p" in trAnnoEnt:
					trAnnoEnt['p3'] = BedAnno.P1toP3(trAnnoEnt['p'])
				if "standard_pHGVS" in trAnnoEnt:
					trAnnoEnt['standard_p3'] = BedAnno.P1toP3(trAnnoEnt['standard_pHGVS'])
				if "alt_pHGVS" in trAnnoEnt:
					trAnnoEnt['alt_p3'] = BedAnno.P1toP3(trAnnoEnt['alt_pHGVS'])

				trAnnoEnt['trVarName'] = BedAnno._getTrVarName(tid, trAnnoEnt)

		annoEnt.var.varName = BedAnno.decide_major(annoEnt)

		return annoEnt

	'''
	=head2 decide_major

		About   : In finalise step, decide a major transcript to report for a var.
		Usage   : my $majorTranscriptVarName = decide_major($annoEnt);
		Returns : A string in the following format:
				  If the transcript has a pName: 
					  mrnaAcc(geneSymbol): cName (pName), 
					  e.g. NM_145651.2(SCGB1C1): c.13C>T (p.R5C)
				  If the transcript does not have a pName: 
					  mrnaAcc(geneSymbol): cName
				  If only intergenic
					  chr: gName (intergenic)
					  e.g. chrX: g.220025A>T (intergenic)
		Note    : use the primaryTag to find reference standard or primary transcript,
				  if only one have the primaryTag "Y", then use it,
				  otherwise, sort by GenePart: 
					1.  CDS
					2.  span
					3.  five_prime_cis_splice_site
					4.  three_prime_cis_splice_site
					5.  ncRNA
					6.  five_prime_UTR
					7.  three_prime_UTR
					8.  interior_intron
					9.  five_prime_UTR_intron
					10. three_prime_UTR_intron
					11. abnormal-intron
					12. promoter
					13. annotation-fail
					14. intergenic_region
				  and choose the first one, if more than one transcript have the 
				  same reference standard and same GenePart, then choose the first
				  one which prior in name sorting.

	=cut
	'''
	@staticmethod
	def decide_major(annoEnt):
		if not hasattr(annoEnt, "trInfo") or 0 == len(annoEnt.trInfo.keys()) or "" in annoEnt.trInfo:
			gHGVS = annoEnt.var.standard_gHGVS if hasattr(annoEnt.var, "standard_gHGVS") else annoEnt.var.gHGVS
			intergenic = "{}{}{}".format(annoEnt.var.chr, ": ", gHGVS)
			if re.search(r'^[\dXYM]', intergenic):
				intergenic = "{}{}".format("chr", intergenic)
			intergenic = "{}{}".format(intergenic, " (intergenic)")
			if re.search(r'^chrMT', intergenic):
				intergenic = re.sub(r'^chrMT', BedAnno.CURRENT_MT, intergenic)
			return intergenic
		else:
			prTrs = dict()
			nonPrTrs = dict()
			for tid in annoEnt.trInfo.keys():
				if "genepart" not in annoEnt.trInfo[tid]:
					raise "Error: cannot decide major transcript before finaliseAnno"
				if str(annoEnt.trInfo[tid]['primaryTag']) == "Y":
					prTrs[tid] = annoEnt.trInfo[tid]['genepart']
				else:
					nonPrTrs[tid] = annoEnt.trInfo[tid]['genepart']

			majorTr = None
			if 0 < len(prTrs.keys()):  # hit primary transcripts
				if 1 == len(prTrs.keys()):
					majorTr = prTrs.keys()[0]
				else:
					majorTr = BedAnno._get_first_tr(prTrs)
			else:  # hit only non-primary transcripts
				majorTr = BedAnno._get_first_tr(nonPrTrs)
			rTrEnt = annoEnt.trInfo[majorTr]

			return annoEnt.trInfo[majorTr]['trVarName']

	# involke parse_annoent to assign detail information to annodb.
	@staticmethod
	def assign_detail(rannodb_k):
		detail = dict()
		for annoblk in sorted(rannodb_k['annos'].keys()):
			offset = rannodb_k['annos'][annoblk]
			tid, ranno = BedAnno.parse_annoent(annoblk)
			detail[tid] = ranno
			detail[tid]['offset'] = offset
		rannodb_k['detail'] = detail
		return rannodb_k

	# just parse annoents when mutation hits.
	@staticmethod
	def parse_annoent(annoent):
		annoinfo = dict()
		# tid         gsym  gid strd blka gpgo exin nsta nsto csta csto wlen mismatch    pr
		# NM_015658.3|NOC2L|26155|-|3U1E|205|EX19E|2817|2801|*508|*492|0|D,879582,879582,.|Y
		infos = re.split(r'\|', annoent)
		if 14 != len(infos):
			raise "".format("Error format of anno ents [", annoent, "]")

		tid = infos.pop(0)
		tags = ["gsym", "gid", "strd", "blka", "gpSO", "exin", "nsta", "nsto", "csta", "csto", "wlen", "mismatch", "pr"]
		for k, v in zip(tags, infos):
			annoinfo[k] = v
		if annoinfo['gpSO'] == 'abnormal_intron':
			annoinfo['gpSO'] = 'abnormal-intron'
		return (tid, annoinfo)


	@staticmethod
	def rev_comp(Seq):
		Seq = Seq.replace('A', '{A}').replace('T', '{T}').replace('C', '{C}').replace('G', '{G}')
		Seq = Seq.format(A='T', T='A', C='G', G='C')[::-1]  # only deal with 'A', 'T', 'G', 'C'
		return Seq

	'''
	=head2 cmpPos

		About   : judge the order for p1 and p2, because the insertion
				  will have a reverted order of left & right position
		Usage   : my $cmpRst = BedAnno->cmpPos($p1, $p2);
		Args    : hgvs positio p1 and p2, with out 'c.' or 'n.' flag
		Return  : 0 for same, 1 for normal order, -1 for reversed order.

	=cut
	'''
	@staticmethod
	def cmpPos(p1, p2):
		if str(p1) == str(p2):
			return 0
		s1, s2, anc1, anc2, int_s1, int_s2, ofst1, ofst2 = None, None, None, None, None, None, None, None

		order_s = {'-u': 1, '-': 2, '': 3, '+': 4, '*': 5, '+d': 6}

		m = re.match(r'([\-\+\*]?)(\d+)([\+\-]?[ud]?)(\d*)$', str(p1))
		if m:
			s1 = m.group(1)
			anc1 = m.group(2)
			int_s1 = m.group(3)
			ofst1 = m.group(4)

		m = re.match(r'([\-\+\*]?)(\d+)([\+\-]?[ud]?)(\d*)$', str(p2))
		if m:
			s2 = m.group(1)
			anc2 = m.group(2)
			int_s2 = m.group(3)
			ofst2 = m.group(4)
		if order_s[s1] < order_s[s2]:
			return 1
		elif order_s[s1] > order_s[s2]:
			return -1
		else:  # $s1 eq $s2
			if str(s1) == '-':
				if int(anc1) > int(anc2):
					return 1
				elif int(anc1) < int(anc2):
					return -1
			else:
				if int(anc1) < int(anc2):
					return 1
				elif int(anc1) > int(anc2):
					return -1
			# anc1 eq anc2
			if order_s[int_s1] < order_s[int_s2]:
				return 1
			elif order_s[int_s1] > order_s[int_s2]:
				return -1
			else:  # int_s2 same
				if re.search(r'\-', str(int_s1)):
					if int(ofst1) > int(ofst2):
						return 1
					elif int(ofst1) < int(ofst2):
						return -1
				else:
					if int(ofst1) < int(ofst2):
						return 1
					elif int(ofst1) > int(ofst2):
						return -1
		return 0


	@staticmethod
	def _getCodingSeq(trdbEnt):
		if "csta" not in trdbEnt:
			return ""
		codingSeq = trdbEnt['seq'][int(trdbEnt['csta']):int(trdbEnt['csto']) ]
		if "cfs" not in trdbEnt:
			return codingSeq
		for fsld in sorted(trdbEnt['cfs'].keys(), reverse=True):
			if (fsld <= 0) or int(fsld) + int(trdbEnt['cfs'][fsld]) <= 0 or fsld >= (len(codingSeq) - 1):
				raise "{}{}{}".format("Error: [", trdbEnt['gene'], "] frameshift position error.")
			if trdbEnt['cfs'][fsld] < 0:
				codingSeq = "{}{}{}".format(codingSeq[:fsld], codingSeq[(int(fsld) + int(trdbEnt['cfs'][fsld])):(
				abs(int(trdbEnt['cfs'][fsld])) + int(trdbEnt['cfs'][fsld]))], codingSeq[fsld:])
			else:
				codingSeq = "{}{}".format(codingSeq[:int(fsld)], codingSeq[(int(fsld) + int(trdbEnt['cfs'][fsld])):])
		return codingSeq

	'''
	=head2 getTrRef

		About   : generate concatenated transcript originated reference
		Usage   : my $trRef = getTrRef( $trannoEnt, $refgenome, $trSeq, $strd );
		Args    : trannoEnt - BedAnno::Anno->{trInfo}->{$tid}
				  refgenome - Unified reference in BedAnno::Var
				  trSeq     - whole transcript
				  strd      - strand of transcript.
		Returns : transcript originated reference
		Notes   : Using sequence from transcript as the exon part,
				  and using the sequence from reference genome
				  as the intron part. and concatenate them.

	=cut
	'''
	@staticmethod
	def getTrRef(trannoEnt, refgenome, trSeq, strd):
		tag_sort = sorted(trannoEnt['trRefComp'].keys(), cmp=BedAnno.trRefSort)
		trRef = ""
		trStart = BedAnno.getTrStart(trannoEnt['rnaBegin'])
		for exin in tag_sort:
			if not re.match(r'EX', exin):
				if refgenome == '=':
					return '='
				int_seq = None
				if 2 == len(trannoEnt['trRefComp'][exin]):
					int_seq = refgenome[int(trannoEnt['trRefComp'][exin][0]):int(trannoEnt['trRefComp'][exin][1])]
				else:
					int_seq = ""
				if not strd:
					int_seq = BedAnno.rev_comp(int_seq)
				trRef = "{}{}".format(trRef, int_seq)
			else:
				trRef = "{}{}".format(trRef , trSeq[int(trStart):(int(trStart) + int(trannoEnt['trRefComp'][exin]))])
				trStart = int(trStart) + int(trannoEnt['trRefComp'][exin])
		return trRef


	@staticmethod
	def trRefSort(a, b):
		atag, anum, btag, bnum = None, None, None, None
		m = re.search(r'^(\D+)(\d+)', a)
		if m:
			atag = m.group(1)
			anum = m.group(2)
		m = re.search(r'^(\D+)(\d+)', b)
		if m:
			btag = m.group(1)
			bnum = m.group(2)
		if (atag is None) or (anum is None) or (btag is None) or (bnum is None):
			raise "{}{}{}{}{}".format("Error exon intron number [", a, ", ", b, "]")
		return cmp(anum, bnum) or cmp(atag, btag)

	# get the next or current on-transcript position
	# from any nDot format HGVS position string.
	@staticmethod
	def getTrStart(trStart):
		m1 = re.search(r'^(\d+)\-', str(trStart))
		m2 = re.search(r'^(\d+)\+', str(trStart))
		m3 = re.search(r'^(\d+)$', str(trStart))
		if re.search(r'^-', str(trStart)):
			return 0
		elif m1:
			return int(m1.group(1)) - 1
		elif m2:
			return int(m2.group(1))
		elif m3:
			return int(m3.group(1)) - 1
		else:
			raise "Error: not available trStart [" + str(trStart) + "]"

	'''
	=head2 _getPairedCodon

		About   : get codon position, codon string, aa string, and frame info
				  for a pair of transcript position
		Usage   : my ($rcinfo5, $rcinfo3) = _getPairedCodon( $trdbEnt, $p5, $p3 );
		Args    : trdbEnt is a sub hash in trInfodb which contains csta, csto
				  for cds start/end position, besides many other feature tags.
				  p5 and p3 is a pair of positions which give a region or single
				  position or maybe an insertion anchor.
		Returns : codon info array ref for p5 and p3. The array ref:
				  [
					  AA position  - 0 for not in cds region.
					  codon string - codon string, e.g. "ATA".
					  aa char      - AA code, 1 bp mode, e.g. "E".
					  polar        - Polar properties, e.g. "P+".
					  frame        - current position's frame info,
									 -1 for not available or in fs-site.
			  [frame-alt]  - frame info around fs-site.
				  ]

	=cut
	'''
	@staticmethod
	def _getPairedCodon(trdbEnt, p5, p3):
		c5 = BedAnno._getCPosFrame(trdbEnt, p5)
		c3 = BedAnno._getCPosFrame(trdbEnt, p3)
		if c5[0] == 1 and c3[0] == 1:  # all normal
			if "nfs" in trdbEnt:
				contain_fs = 0
				for fsld in trdbEnt['nfs'].keys():
					if p5 <= fsld and p3 >= fsld + trdbEnt['nfs'][fsld]:
						contain_fs = 1
						break
				if contain_fs == 1:
					return ([c5[1][0]]+BedAnno._genCodonInfo(trdbEnt, c5[1][0])+[ -1, c5[1][1]],
							[c3[1][0]]+BedAnno._genCodonInfo(trdbEnt, c3[1][0])+[ -1, c3[1][1]])
			return ([c5[1][0]]+BedAnno._genCodonInfo(trdbEnt, c5[1][0])+ [c5[1][1]],
					[c3[1][0]]+BedAnno._genCodonInfo(trdbEnt, c3[1][0])+[ c3[1][1]])
		elif c5[0] == 0 and c3[0] == 0:  # all in deleted frame
			if c5[1][0] == c3[1][0]:  # same frame
				return ([c5[2][0]]+ BedAnno._genCodonInfo(trdbEnt, c5[2][0])+[ -1],
						[c5[1][0]]+ BedAnno._genCodonInfo(trdbEnt, c5[1][0])+[ -1])
			else:  # different frames
				return ([c5[2][0]]+ BedAnno._genCodonInfo(trdbEnt, c5[2][0])+[ -1, c5[2][1]],
						[c3[1][0]]+ BedAnno._genCodonInfo(trdbEnt, c3[1][0])+[ -1, c3[1][1]])
		elif c5[0] == 2 and c3[0] == 2:  # all in duplicated frame
			if BedAnno._chkCodonPosRel(c5, c3):   # connected
				return ([c5[1][0]]+ BedAnno._genCodonInfo(trdbEnt, c5[1][0])+[ -1, c5[1][1]],
					[c3[2][0]]+ BedAnno._genCodonInfo(trdbEnt, c3[2][0])+[ -1, c3[2][1]])
			else:  # non-connected, this case don't exists in current db,
				# we will only take the first two returned value
				# and give up the latter two.
				return ([c5[1][0]]+BedAnno._genCodonInfo(trdbEnt, c5[1][0])+[ -1, c5[1][1]],
						[c3[1][0]]+BedAnno._genCodonInfo(trdbEnt, c3[1][0])+[ -1, c3[1][1]],
						[c5[2][0]]+BedAnno._genCodonInfo(trdbEnt, c5[2][0])+[ -1, c5[2][1]],
						[c3[2][0]]+ BedAnno._genCodonInfo(trdbEnt, c3[2][0])+[ -1, c3[2][1]])

		else:
			r5_ret, r3_ret = None, None
			if c5[0] == 1:
				r5_ret = [c5[1][0]]+ BedAnno._genCodonInfo(trdbEnt, c5[1][0])+ [c5[1][1]]
			elif c5[0] < 0:
				r5_ret = [0, "", "", "", -1]
			elif c5[0] == 0:
				r5_ret = [c5[2][0]]+BedAnno._genCodonInfo(trdbEnt, c5[2][0])+[ -1, c5[2][1]]
			elif c5[0] == 2:
				r5_ret = [c5[1][0]]+ BedAnno._genCodonInfo(trdbEnt, c5[1][0])+[ -1, c5[1][1]]
			if c3[0] == 1:
				r3_ret = [c3[-1][0]]+ BedAnno._genCodonInfo(trdbEnt, c3[-1][0])+[ c3[-1][1]]
			elif c3[0] < 0:
				r3_ret = [0, "", "", "", -1]
			elif c3[0] == 0:
				r3_ret = [c3[1][0]]+ BedAnno._genCodonInfo(trdbEnt, c3[1][0])+[ c3[1][1]]
			elif c3[0] == 2:
				r3_ret = [c3[2][0]]+ BedAnno._genCodonInfo(trdbEnt, c3[2][0])+[ c3[2][1]]
			return (r5_ret, r3_ret)


	@staticmethod
	def _getCPosFrame(trdbEnt, p):
		if "csta" not in trdbEnt or not re.search(r'^\d+$', str(p)):
			return [-1]
		cds_p = int(p) - int(trdbEnt['csta'])
		return BedAnno._getCPosFrame_by_cdsPos(trdbEnt, cds_p)

	'''
	=head2 _getCPosFrame_by_cdsPos

		About   : calculate codon position and frame info by
				  involving frameshift case in consideration.
		Usage   : my @cPosFrame = _getCPosFrame_by_cdsPos( $trdbEnt, $cds_p );
		Args    : cds_p is relative position to the first bp in start codon on trans
		Returns : an array in the format of:
				  ( $stat, $ra_posrefs1, $ra_posrefs2 )
				  ra_posrefs* is an array ref for the following array
				  [ $codonPos, $frame ]

				  for non-cds region case stat is -1, without ra_posrefs.
				  for normal case stat is 1, 
					  ra_posrefs1 is the information of input position,
					  without ra_posrefs2 specified.
				  for deleted frame case stat is 0,
					  ra_posrefs1 is the information of position on 5' 
					  side next to the deleted frame, and ra_posrefs2
					  is on 3' side.
				  for duplicated frame case stat is 2,
					  ra_posrefs1 is the 5' most hit of input position,
					  ra_posrefs2 is the corresponding position on the 
					  3' duplicated region.

	=cut
	'''
	# give codon number and frame with involving frameshift case
	@staticmethod
	def _getCPosFrame_by_cdsPos(trdbEnt, p):
		p = int(p)
		if "csta" not in trdbEnt or trdbEnt['csta'] == '.' or not re.search(r'^\d+$', str(p)) or int(p) <= 0 or int(p) > (int(trdbEnt['csto']) - int(trdbEnt['csta'])):
			return [-1]
		if "cfs" in trdbEnt:
			total_fs = 0
			for fsld in sorted(trdbEnt['cfs'].keys(), cmp=lambda a, b: cmp(a, b)):
				if int(fsld) < p and int(fsld) + int(trdbEnt['cfs'][fsld]) < p:
					# normal
					total_fs += int(trdbEnt['cfs'][fsld])
				elif int(fsld) < p and int(fsld) + int(trdbEnt['cfs'][fsld]) >= p:
					# in deleted frames
					codon_pre = fsld - int(total_fs)
					pP_pre, frame_pre = BedAnno._calPosFrame(codon_pre)
					pP_lat, frame_lat = None, None
					if frame_pre < 2:
						pP_lat = pP_pre
						frame_lat = int(frame_pre) + 1
					else:
						pP_lat = int(pP_pre) + 1
						frame_lat = 0
					return [0, [pP_pre, frame_pre], [pP_lat, frame_lat]]
				elif int(fsld) >= p and int(fsld) + int(trdbEnt['cfs'][fsld]) < p:
					# in dup frames
					codon_1 = p - total_fs
					codon_2 = p - total_fs + abs(int(trdbEnt['cfs'][fsld]))
					pP_1, frame_1 = BedAnno._calPosFrame(codon_1)
					pP_2, frame_2 = BedAnno._calPosFrame(codon_2)
					return [2, [pP_1, frame_1], [pP_2, frame_2]]
				else:  # after the current position
					break
			p -= total_fs
		return [1, BedAnno._calPosFrame(p)]

	# calculate codon number and frame directly be cds position
	@staticmethod
	def _calPosFrame(cds_p):
		pP = int(cds_p) / 3
		if int(cds_p) % 3 > 0:
			pP += 1
		frame = 2 - (pP * 3 - int(cds_p))
		return [pP, frame]

	'''
	=head2 trWalker

		About   : walk around the variant position to find possible
				  repeat start/end, and return the recalculated
				  trBegin and trEnd, together with the real
				  transcript originated variants and unified property
				  Current implementation won't walk around in the following
				  cases:
				  1. no-call
				  2. annotation-fail
				  3. snv or mnp
				  4. non-exon region
				  5. span case or any edge case
				  6. delins without repeat.
		Usage   : ($trBegin, $trEnd, $real_var) = 
				  $beda->trWalker($tid, $rtrinfo);
		Args    : trAcc and hash ref of trInfo annotation for trAcc,
			  

	=cut
	'''
	def trWalker(self, tid, rtrinfo):
		# reparse the transcript originated var
		real_var = BedAnnoVar.BedAnnoVar(tid, 0, len(rtrinfo['trRef']), rtrinfo['trRef'], rtrinfo['trAlt'])
		Unified = real_var.getUnifiedVar('+')
		real_p, real_r, real_a, real_rl, real_al = Unified

		trBegin = BedAnno.reCalTrPos_by_ofst(rtrinfo, real_p)
		trEnd = BedAnno.reCalTrPos_by_ofst(rtrinfo, (real_p + real_rl - 1))

		# real_rl == real_al    skip snv and mnp
		# rtrinfo['ei_Begin']!=rtrinfo['ei_End']    skip span case or any edge case
		# not re.search(r'^\d+$',trBegin) or not re.search(r'^\d+$',trEnd)     skip non exon begin/end position
		# real_r.find(real_a)!=0 and real_a.find(real_r)!=0    skip delins without repeat
		if (real_rl == real_al) or (rtrinfo['ei_Begin'] != rtrinfo['ei_End']) or (
			not re.search(r'^\d+$', str(trBegin)) or not re.search(r'^\d+$', str(trEnd))) or (
				real_r.find(real_a) != 0 and real_a.find(real_r) != 0):
			# no correction
			return (trBegin, trEnd, real_var, Unified)

		qtid = tid
		qtid = re.sub(r'\-\d+$', '', str(qtid))
		trSeq = self.trInfodb[qtid]['seq']
		ref_sta, ref_sto, trRef, trAlt = BedAnno.walker(trBegin, trEnd, trSeq, real_r, real_a, real_rl, real_al)

		if int(ref_sta) != int(trBegin) or int(ref_sto) != int(trEnd):
			trBegin = ref_sta
			trEnd = ref_sto
			real_var = BedAnnoVar.BedAnnoVar(tid, 0, len(trRef), trRef, trAlt)
			Unified = real_var.getUnifiedVar('+')

		return (trBegin, trEnd, real_var, Unified)

	# new var generated from new() in BedAnnoVar the real transcript var
	# will bring some new transcript start/end information
	# here will recalculate it from the offset.
	@staticmethod
	def reCalTrPos_by_ofst(trannoEnt, trRef_ofst):
		if int(trRef_ofst) == 0:
			return trannoEnt['rnaBegin']
		if int(trRef_ofst) == -1 and "preStart" in trannoEnt:
			return trannoEnt['preStart']['nDot']

		tag_sort = sorted(trannoEnt['trRefComp'].keys())
		cumulate_len = 0
		cur_ex_start, intOfst = None, None
		m = re.search(r'^(\-?\d+)([\+\-]?\d*)$', str(trannoEnt['rnaBegin']))
		if m:
			cur_ex_start = m.group(1)
			intOfst = m.group(2)
		if re.search(r'^\+', str(intOfst)):
			cur_ex_start += 1
		if re.search(r'^-', str(cur_ex_start)):
			intOfst = cur_ex_start
			cur_ex_start = 1
		if str(intOfst) == '':
			intOfst = 0

		for exin in tag_sort:
			cur_blk_len = None
			if not re.search(r'^EX', str(exin)):
				if 2 > len(trannoEnt['trRefComp'][exin]):
					cur_blk_len = 0
				else:
					cur_blk_len = int(trannoEnt['trRefComp'][exin][1]) - int(trannoEnt['trRefComp'][exin][0])
			else:
				cur_blk_len = trannoEnt['trRefComp'][exin]
			cur_blk_ofst = (trRef_ofst - cumulate_len)
			if int(cur_blk_ofst) >= int(cur_blk_len):
				if int(cur_blk_ofst) == int(cur_blk_len) and str(exin) == str(tag_sort[-1]):
					return trannoEnt['postEnd']['nDot']
				else:
					cumulate_len += int(cur_blk_len)
					if re.search(r'^EX', exin):  # only cumulate ex pos
						cur_ex_start =int(cur_ex_start)+ int(cur_blk_len)
						intOfst = 0
			else:  # hit current block
				if not re.search(r'^EX', str(exin)):
					if int(intOfst) < 0 and int(cur_ex_start) == 1:  # in promoter
						if int(cumulate_len) > 0:
							raise "{}{}{}".format("Error: non-zero cumulate_len in promoter [", cumulate_len, "].")
						return (int(intOfst) + int(cur_blk_ofst))
					elif re.search(r'^Z', str(exin)):
						return "{}{}".format('+', int(cur_blk_ofst) + 1)
					else:  # intron
						blk_begin, blk_end = None, None
						if str(exin) == str(tag_sort[0]):
							blk_begin = trannoEnt['rnaBegin']
						else:
							blk_begin = "{}{}".format(cur_ex_start - 1, "+1")
						if str(exin) == str(tag_sort[-1]):
							blk_end = trannoEnt['rnaEnd']
						else:
							blk_end = "{}{}".format(cur_ex_start, "-1")
						return BedAnno.getIntrPos(blk_begin, blk_end, cur_blk_ofst, cur_blk_len)
				else:
					return (int(cur_ex_start) + int(cur_blk_ofst))
		raise "Error: out of trRefComp range [" + str(trRef_ofst) + "]."


	@staticmethod
	def walker(ref_sta, ref_sto, whole_seq, ref, alt, reflen, altlen):
		seqlen = len(whole_seq)

		ori_walker = None
		track_opt = None

		if reflen < altlen:
			ref_sta = ref_sto + 1
			ori_walker = alt[reflen:]
			track_opt = 0  # walk on alt
		else:
			ref_sta += altlen
			ori_walker = ref[altlen:]
			track_opt = 1  # walk on ref

		# walk to 3' most
		cur_walker = list(ori_walker)
		walk_forward_step = 0
		for p in range(ref_sto, seqlen):
			if whole_seq[p:p + 1] == cur_walker[0]:
				cur_walker.append(cur_walker.pop(0))
				if p == seqlen - 1:
					walk_forward_step = seqlen - ref_sto
					break
			else:
				walk_forward_step = p - ref_sto
				break

		ref_sto += walk_forward_step  # the final end on reference track

		# use the 3' most element as the final matched difference
		match_target = "".join(cur_walker)

		forward_footprint = whole_seq[(ref_sta - 1):ref_sta - 1 + walk_forward_step]

		# walk to 5' most
		cur_walker = list(ori_walker)
		walk_back_step = 0
		for q in range(ref_sta - 1, 0, -1):
			if whole_seq[(q - 1):q] == cur_walker[-1]:
				cur_walker = list(cur_walker[-1]) + cur_walker[:-1]
				if q == 1:
					walk_back_step = ref_sta - 1
					break
			else:
				walk_back_step = ref_sta - q - 1
				break

		back_most = ref_sta - walk_back_step
		backward_footprint = whole_seq[(ref_sta - walk_back_step - 1):(ref_sta - 1)]

		cur_short_track = "{}{}".format(backward_footprint, forward_footprint)
		cur_long_track = "{}{}".format(cur_short_track, match_target)
		target_sta = cur_long_track.index(match_target)
		ref_sta = back_most + target_sta

		renew_ref, renew_alt = None, None
		if track_opt:
			renew_ref = cur_long_track[target_sta:]
			renew_alt = cur_short_track[target_sta:]
		else:
			renew_ref = cur_short_track[target_sta:]
			renew_alt = cur_long_track[target_sta:]

		return (ref_sta, ref_sto, renew_ref, renew_alt)

	# calculate position in intron region
	@staticmethod
	def getIntrPos(begin, end, ofst, blk_len):
		exPosBegin, exPosEnd, intPosBegin, intPosEnd = None, None, None, None
		m = re.search(r'^(\d+)([\+\-]\d+)$', str(begin))
		if m:
			exPosBegin = m.group(1)
			intPosBegin = m.group(2)
		m = re.search(r'^(\d+)([\+\-]\d+)$', str(end))
		if m:
			exPosEnd = m.group(1)
			intPosEnd = m.group(2)
		if int(exPosBegin) == int(exPosEnd):
			if int(intPosBegin) < 0:
				return "{}{}".format(exPosBegin,  int(intPosBegin) + int(ofst))
			else:
				return "{}{}{}".format(exPosBegin, '+', int(intPosBegin) + int(ofst))
		else:
			wlen = int(intPosBegin) - int(intPosEnd) + int(blk_len) - 2
			lofst = int(intPosBegin) + int(ofst) - 1
			if lofst < (wlen / 2.0 - 1):
				return ("{}{}{}".format(exPosBegin, '+', lofst + 1))
			else:
				return ("".format(exPosEnd, '-', wlen - lofst))

	'''
	=head2 P1toP3

		About : Change 1 letter format of pHGVS string to 3 letter format
		Usage : my $p3 = P1toP3($p1);

	=cut
	'''
	@staticmethod
	def P1toP3(p1):
		m1 = re.search(r'^p\.([A-Z\*])(\d+)([A-Z\*])((fs\*.+)?)$', p1)
		m2 = re.search(r'^p\.([A-Z\*])(\d+)(del|dup|\[.*\])$', p1)
		m3 = re.search(r'^p\.([A-Z\*])(\d+)delins([A-Z\*]+)$', p1)
		m4 = re.search(r'^p\.([A-Z\*])(\d+)_([A-Z\*])(\d+)(del|dup|\[.*\])$', p1)
		m5 = re.search(r'^p\.([A-Z\*])(\d+)_([A-Z\*])(\d+)((del)?ins)([A-Z\*]+)$', p1)
		m6 = re.search(r'^p\.([A-Z\*])(\d+)((delins)?)\?$', p1)
		m7 = re.search(r'^p\.([A-Z\*])(\d+)_([A-Z\*])(\d+)(delins\?)$', p1)
		m8 = re.search(r'^p\.\(([A-Z\*])(\d+)=\)$', p1)
		if m1:
			# missense, frameshift
			return "{}{}{}{}{}".format('p.', BedAnno.C1toC3[m1.group(1)], m1.group(2), BedAnno.C1toC3[m1.group(3)],
									   m1.group(4))
		elif m2:
			# 1 bp deletion, duplication, repeats
			return "{}{}{}{}".format('p.', BedAnno.C1toC3[m2.group(1)], m2.group(2), m2.group(3))
		elif m3:
			# 1 bp delins
			former = "{}{}{}{}".format('p.', BedAnno.C1toC3[m3.group(1)], m3.group(2), 'delins')
			singles = list(m3.group(3))
			latter = "".join([BedAnno.C1toC3[k] for k in singles])
			return "{}{}".format(former, latter)
		elif m4:
			# long deletion, duplication, repeats
			return "{}{}{}{}{}{}{}".format('p.', BedAnno.C1toC3[m4.group(1)], m4.group(2), '_', BedAnno.C1toC3[m4.group(3)],
										   m4.group(4), m4.group(5))
		elif m5:
			# insertion, long delins
			former = "{}{}{}{}{}{}{}".format('p.', BedAnno.C1toC3[m5.group(1)], m5.group(2), '_',
											 BedAnno.C1toC3[m5.group(3)], m5.group(4), m5.group(5))
			singles = list(m5.group(7))
			latter = "".join([BedAnno.C1toC3[k] for k in singles])
			return "{}{}".format(former, latter)
		elif m6:
			# 1 base no-call
			return "{}{}{}{}{}".format('p.', BedAnno.C1toC3[m6.group(1)], m6.group(2), m6.group(3), '?')
		elif m7:
			# long no-call
			return "{}{}{}{}{}{}{}".format('p.', BedAnno.C1toC3[m7.group(1)], m7.group(2), '_', BedAnno.C1toC3[m7.group(3)],
										   m7.group(4), m7.group(5))
		elif m8:
			return "{}{}{}{}".format('p.(', BedAnno.C1toC3[m8.group(1)], m8.group(2), '=)')
		else:
			return p1


	@staticmethod
	def _checkDistEdge(tranno, dbk):
		edge_pos = None
		m = re.search(r'^(\d+)[\+\-]', dbk['nsta'])
		if m:
			edge_pos = m.group(1)
		if edge_pos is not None:
			if (re.search(r'^\d+$', str(tranno['rnaBegin'])) and abs(tranno['rnaBegin'] - int(edge_pos)) < 3) or (
				re.search(r'^\d+$', str(tranno['rnaEnd'])) and abs(tranno['rnaEnd'] - int(edge_pos)) < 3):
				return 1
		return 0


	@staticmethod
	def _getTrVarName(trID, rTrEnt):
		trhit = None

		if re.search(r'N[MR]_MT-', trID):
			trhit = rTrEnt['geneSym']
		else:
			trhit = "{}{}{}{}".format(trID, "(", rTrEnt['geneSym'], ")")

		if rTrEnt['func'] == 'annotation-fail':
			trhit = "{}{}".format(trhit, ": annotation-fail")
		else:
			cHGVS = rTrEnt['standard_cHGVS'] if "standard_cHGVS" in rTrEnt else rTrEnt['c']
			trhit = "{}{}{}".format(trhit, ": ", cHGVS)

		if "p" in rTrEnt and rTrEnt['p'] != "":
			pHGVS = rTrEnt['standard_pHGVS'] if "standard_pHGVS" in rTrEnt else rTrEnt['p']
			trhit = "{}{}{}{}".format(trhit, " (", pHGVS, ")")
		return trhit


	@staticmethod
	def _get_first_tr(rGeneParts_hash):
		def sort_f(a, b):
			if cmp(BedAnno.GenePartsOrder[rGeneParts_hash[a]], BedAnno.GenePartsOrder[rGeneParts_hash[b]]) == -1:
				return -1
			elif cmp(BedAnno.GenePartsOrder[rGeneParts_hash[a]], BedAnno.GenePartsOrder[rGeneParts_hash[b]]) == 1:
				return 1
			else:
				return cmp(a, b)

		ordered_Trs = sorted(rGeneParts_hash.keys(), cmp=sort_f)
		return ordered_Trs[0]

	# The standard gHGVS and pHGVS can be generated directly from
	# original HGVS string, but cHGVS need more specification
	@staticmethod
	def gen_standard_gphgvs(rep_hgvs):
		m = re.search(r'^([gm]\.)(\d+)([ACGTN]+)\[(\d+)>(\d+)\]$', rep_hgvs)
		m2 = re.search(r'^p\.(\D)(\d+)((_\D)(\d+))?\[(\d+)>(\d+)\]$', rep_hgvs)
		if m:
			sym, anchor, rep, ref_cn, alt_cn = m.group(1), m.group(2), m.group(3), m.group(4), m.group(5)
			replen = len(rep)
			if int(alt_cn) - int(ref_cn) == 1:  # convert to dup
				if int(replen) == 1:
					return "{}{}{}{}".format(sym, int(anchor) + int(ref_cn) - 1, 'dup', rep)
				else:
					return "{}{}{}{}{}{}".format(sym, int(anchor) + int(ref_cn - 1) * int(replen), '_', int(anchor) + int(ref_cn) * int(replen) - 1,
												 'dup', rep)
			else:
				return "{}{}{}{}{}{}".format(sym, anchor, rep, '[', alt_cn, ']')
		elif m2:
			anchor1_char, anchor1 = m2.group(1), m2.group(2)
			anchor2_char, anchor2 = m2.group(4), m2.group(5)
			ref_cn, alt_cn = m2.group(6), m2.group(7)
			if not anchor2_char:
				anchor2_char = ""
			if not anchor2:
				anchor2 = ""

			if int(alt_cn) - int(ref_cn) == 1:
				if not m2.group(3) or str(m2.group(3)) == '':
					return "{}{}{}{}".format('p.', anchor1_char, int(anchor1) + int(ref_cn) - 1, 'dup')
				else:
					replen = int(anchor2) - int(anchor1) + 1
					return "{}{}{}{}{}{}".format('p.', anchor1_char, int(anchor1) + (int(ref_cn) - 1) * int(replen), anchor2_char,
												 int(anchor1) + int(ref_cn) * int(replen) - 1, 'dup')
			else:
				return "{}{}{}{}{}{}{}{}".format("p.", anchor1_char, anchor1, anchor2_char, anchor2, '[', alt_cn, ']')
		else:
			return rep_hgvs

	# only ghgvs can generate alt HGVS directly from repeat format string
	@staticmethod
	def gen_alt_ghgvs(rep_hgvs):
		m = re.search(r'^([gm]\.)(\d+)([ACGTN]+)\[(\d+)>(\d+)\]$', rep_hgvs)
		if m:
			sym, anchor, rep, ref_cn, alt_cn = m.group(1), int(m.group(2)), m.group(3), int(m.group(4)), int(m.group(5))
			replen = len(rep)

			if ref_cn > alt_cn:  # deletion
				if ref_cn - alt_cn == 1 and replen == 1:
					return "{}{}{}{}".format(sym, anchor + ref_cn * replen - 1, 'del', rep)
				else:
					return "{}{}{}{}{}{}".format(sym, anchor + alt_cn * replen, '_', anchor + ref_cn * replen - 1, 'del',
												 rep * (ref_cn - alt_cn))
			else:  # insertion
				return "{}{}{}{}{}{}".format(sym, anchor + ref_cn * replen - 1, '_', anchor + ref_cn * replen, 'ins',
											 rep * (alt_cn - ref_cn))
		return rep_hgvs

	# let unavailable pr pos to be 0
	@staticmethod
	def _getCoveredProd(trdbEnt, chgvs_5, chgvs_3):
		if re.search(r'^\-', str(chgvs_3)) or re.search(r'^\*',str(chgvs_5)):  # re.search(r'^\-',chgvs_3) 5'UTR re.search(r'^\*',chgvs_5) 3'UTR
			return (0, 0)

		prb, pre = None, None
		prb_anchor, prb_sig = None, None
		pre_anchor, pre_sig = None, None

		m = re.search(r'^(\d+)([\+\-])', str(chgvs_5))
		if m:
			prb_anchor = m.group(1)
			prb_sig = m.group(2)

		m = re.search(r'^(\d+)([\+\-])', str(chgvs_3))
		if m:
			pre_anchor = m.group(1)
			pre_sig = m.group(2)

		if (prb_anchor is not None) and (prb_sig is not None) and (pre_anchor is not None) and (pre_sig is not None):
			if prb_anchor == pre_anchor or ((prb_anchor == int(pre_anchor) - 1) and str(prb_sig) == '+' and str(pre_sig) == '-'):  # prb_anchor == pre_anchor -> assume no 1 bp exon
				return (0, 0)  # single intron

		m = re.search(r'^(\d+)', str(chgvs_5))
		if re.search(r'^\-', str(chgvs_5)):
			prb = 1
		elif m:
			tmp_cdsp = m.group(1)
			prb, tmp_cdsp = BedAnno._calPosFrame(tmp_cdsp)
		else:
			prb = 0

		m = re.search(r'^(\d+)', str(chgvs_3))
		if re.search(r'^\*', str(chgvs_3)):
			pre = int(trdbEnt['plen']) + 1
		elif m:
			tmp_cdsp2 = m.group(1)
			pre, tmp_cdsp2 =BedAnno. _calPosFrame(tmp_cdsp2)
		else:
			pre = 0

		return (prb, pre)

	# generate codon together with AA and Polar by codon number
	@staticmethod
	def _genCodonInfo(trdbEnt, pP):
		if "cseq" not in trdbEnt:
			trdbEnt['cseq'] = BedAnno._getCodingSeq(trdbEnt)
		if int(pP) <= 0 or int(pP) - 1 > int(trdbEnt['plen']):
			raise "".format("Error: [", trdbEnt['gene'], "] no codon info for ", pP, " position.")
		codon = trdbEnt['cseq'][(pP - 1) * 3:(pP - 1) * 3 + 3]
		rtrans_opts = dict()

		if re.search(r'^MT\-', str(trdbEnt['gene'])):  # chrM gene
			rtrans_opts['mito'] = 1
		if int(pP) <= int(trdbEnt['plen']):  # plen not involve terminal codon
			rtrans_opts['nostop'] = 1
		elif int(pP) == int(trdbEnt['plen']) + 1:  # terminal with polyA complement
			if "A" in trdbEnt:
				codon = "{}{}".format(codon, 'A' * (3 - len(codon)))
		return [codon]+ BedAnno._getAAandPolar(codon, rtrans_opts)


	@staticmethod
	def _getAAandPolar(codon, rtrans_opts):
		aa, zero = BedAnno.translate(codon, rtrans_opts)
		polar = None
		if aa in BedAnno.Polar:
			polar = BedAnno.Polar[aa]
		else:
			if re.search(r'N', str(codon)):
				polar = '?'
			else:
				polar = '.'
		return [aa, polar]

	'''
	=head2 translate

		About   : Translate nucleotides to aa seqs
		Usage   : my ($aa_seq, $next_frame) = translate( $nucl, { mito => 1, polyA => 1 } );
		Args    : first args should be the nucleotide seq to be translated,
				  the second args is a hash ref of optional args (all boolean):
				  mito   : indicate it is for mDNA (mitochondrion)
				  nostop : indicate there's no stop codon in aa seq,
						   translate 'UGA' to 'U', and other stop codon to 'X'
				  polyA  : indicate extra A should be added to 3'end of codon,
						   to help encode a stop codon (usually used with mito).
		Returns : $aa_seq is the aa sequence, 
				  $next_frame gives the next base's frame to the 3'end of sequence.

	=cut
	'''
	@staticmethod
	def translate(nucl, ropt=None):
		mito, nostop, polyA = None, None, None
		if ropt is not None and type(ropt) is dict:

			if "mito" in ropt:
				mito = 1
			else:
				mito = 0
			if "nostop" in ropt:
				nostop = 1
			else:
				nostop = 0
			if "polyA" in ropt:
				polyA = 1
			else:
				polyA = 0

		lenN = len(nucl)
		frame_next = int(lenN) % 3
		if int(frame_next) != 0 and polyA:
			nucl = "{}{}".format(nucl, 'A' * (3 - frame_next))
		nucl = nucl.upper()
		prot = ""
		for i in range(0, len(nucl), 3):
			aa = None
			codon = nucl[i:i + 3]
			if codon in BedAnno.C1:
				aa = BedAnno.C1[codon]
			else:
				if re.search(r'N', str(codon)):
					aa = '?'
				else:
					aa = '.'
			if mito and str(codon) == 'ATA':
				aa = 'M'
			if str(aa) == '*' and mito and str(codon) == 'TGA':
				aa = 'W'
			elif str(aa) == '*' and nostop:
				if str(codon) == 'TGA':
					aa = 'U'
				else:
					aa = 'X'
			elif str(aa) == '*':
				prot = "{}{}".format(prot, aa)
				return (prot, 0)
			prot = "{}{}".format(prot, aa)
		if nostop:
			prot = re.sub(r'[UX]$', '*', str(prot))
		if frame_next != 0:
			prot = re.sub(r'\.$', '', str(prot))
		return (prot, frame_next)

	# change the 1based nDot position into cDot format
	@staticmethod
	def _cPosMark(trPos, csta, csto, trlen):  # 1based	0-based 1based
		if trPos == "":
			return trPos

		cDot = None
		m = re.search(r'^([\+\-]?)(\d+)([\+\-]?)(\d*)$', str(trPos))
		if m:
			s1, anchor, s2, ofst = m.group(1), m.group(2), m.group(3), m.group(4)
			if s1 == '' and s2 == '':
				cDot = int(anchor) - int(csta)
				if cDot <= 0:
					cDot -= 1
				if int(anchor) > int(csto):
					cDot = "{}{}".format('*', int(anchor) - int(csto))
			elif s1 != '' and s2 != '':
				raise "{}{}{}".format("Error: not a valid nDot string: [", trPos, "]")
			elif s1 == '-':
				tmp = None
				if int(csta) == 0:
					tmp = 1
				else:
					tmp = -int(csta)
				cDot = "{}{}{}".format(tmp, '-u', anchor)
			elif s1 == '+':
				tmp = None
				if str(csto) == str(trlen):
					tmp = int(csto) - int(csta)
				else:
					tmp = "{}{}".format('*', int(trlen) - int(csto))
				cDot = "{}{}{}".format(tmp, '+d', anchor)
			elif re.search(r'[\+\-]', s2):
				cDot_tmp = BedAnno._cPosMark(anchor, csta, csto, trlen)
				cDot = "{}{}{}".format(cDot_tmp, s2, ofst)
			else:
				raise "{}{}{}".format("Error: unmatched nDot string: [", trPos, "]")
		return cDot

	@staticmethod
	def pr_r_comp_pr_a(pr_r, pr_a):
		if len(pr_r) >= len(pr_a):
			larger = pr_r
			smaller = pr_a
		else:
			larger = pr_a
			smaller = pr_r
		try:
			result = larger.index(smaller)
		except:
			return -1
		return result

	'''
	=head2 prWalker

		About   : act as trWalker, but on protein sequence
		Usage   : $prVar = $beda->prWalker( $prSeq,  $prBegin, $prEnd, $prRef, $prAlt );
		Args    : prSeq   - whole protein sequence
			  prBegin - protein variant Begin position (1 based)
			  prEnd   - protein variant End position (1 based)
			  prRef   - protein variant reference sequence
			  prAlt   - protein variant alt sequence

	=cut
	'''
	def prWalker(self, prSeq=None, prBegin=None, prEnd=None, prRef=None, prAlt=None):
		prVar = BedAnnoVar.BedAnnoVar("nouse", prBegin - 1, prEnd, prRef, prAlt)
		pr_P, pr_r, pr_a, pr_rl, pr_al = prVar.getUnifiedVar('+')
		if (pr_al is not None and pr_rl == pr_al) or (0 != BedAnno.pr_r_comp_pr_a(pr_r,pr_a)):
			return prVar

		prBegin = pr_P + 1
		prEnd = pr_P + pr_rl

		ref_sta, ref_sto, renew_prRef, renew_prAlt = BedAnno.walker(prBegin, prEnd, prSeq, pr_r, pr_a, pr_rl, pr_al)

		if prBegin != ref_sta or prEnd != ref_sto:
			prVar = BedAnnoVar.BedAnnoVar("nouse", ref_sta - 1, ref_sto, renew_prRef, renew_prAlt)

		return prVar

	'''
	=head2 _cdsubstr

		About   : substr from transcript seqeunce involving frameshift changing.
		Usage   : my $codonStr = _cdsubstr( $trdbEnt, $start, $length, $replace );
		Args    : trdbEnt is a sub hash in trInfodb which contains csta, csto.
				  start is 0 based start position of sub string on transcript seq.
				  length is the length of target region on transcript seq.
				  replace is the alternative.
		Returns : A substring cut from the transcript sequence, with frameshift involved.
				  Behave like function "substr", but for replace mode, 
				  it returns the whole transcript seq after replacement.

	=cut
	'''
	@staticmethod
	def _cdsubstr(trdbEnt=None, start=None, len=None, replace=None):
		if len is None:
			len = int(trdbEnt['len']) - int(start)
		if "csta" not in trdbEnt or "nfs" not in trdbEnt:
			if replace is not None:
				trSeq = trdbEnt['seq']
				trSeq = "{}{}{}".format(trSeq[:start] , replace , trSeq[(start + len):])
				return trSeq
			else:
				return trdbEnt['seq'][start:(start + len)]
		else:
			if start < 0:
				start = trdbEnt['len'] + start
			if len < 0:
				len = trdbEnt['len'] + len - start
			start_ofst = BedAnno._calfsOfst(trdbEnt, start)
			stop_ofst = BedAnno._calfsOfst(trdbEnt, (start + len - 1))

			if "cseq" not in trdbEnt:
				trdbEnt['cseq'] = BedAnno._getCodingSeq(trdbEnt)
			new_wholeseq = "{}{}{}".format(trdbEnt['cseq'][:int(trdbEnt['csta'])], trdbEnt['cseq'],
										   trdbEnt['cseq'][int(trdbEnt['csto']):])

			returned_seq = None
			max_len = -1
			for sta_ofst in start_ofst:
				new_start = start - sta_ofst
				for sto_ofst in stop_ofst:
					new_stop = start + len - 1 - sto_ofst
					if new_start > new_stop:
						continue
					new_stop += 1
					new_len = new_stop - new_start
					if new_len > max_len:
						if replace is not None:
							returned_seq = new_wholeseq
							returned_seq = returned_seq[:new_start] + replace + returned_seq[(new_start + new_len):]
							max_len = new_len
						else:
							returned_seq = new_wholeseq[new_start:(new_start + new_len)]
							max_len = new_len
			# only return the longest substring
			return returned_seq

	'''
	=head2 anno

		About   : Annotate single short variation by annotation db.
		Usage   : my $anno_ent = $beda->anno( 'chr20', 1234567, 1234569, 'AG', 'TGGG' );
				  or $anno_ent = $beda->anno( 'chr20', 1234568, 'AG', 'AGGG' );
		Args    : for CG's shell variants, need 5 args in UCSC coordinates
				  (0-based start), they are:
					chr id, chr start, chr end, reference, alternative.
				  for variants in VCF, need 4 args, which is lack of 
					chr end, and "chr start" is in 1-based coordinates.
				  for crawler: a input object with keys: 
			chr,begin,referenceSequence,variantSequence,[end].
					if end is specified, then use 0-based coordinates,
					otherwise 1-based (VCF) coordinates will be used.
		Returns : a hash ref of annotation informations, see varanno().

	=cut
	'''
	def anno(self, *args):
		debugOpt = None
		if hasattr(self, "debug"):
			debugOpt = 1
		else:
			debugOpt = 0
		t0 = None
		if debugOpt:
			t0 = [time.time()]
		var = BedAnnoVar.BedAnnoVar(*args)
		##var = eval("BedAnnoVar.BedAnnoVar("+",".join(map(str,args))+")")

		t1 = None
		if debugOpt:
			t1 = [time.time()]
			sys.stderr.write("{}{}{}".format("BedAnno->anno [new BedAnno::Var] ... ", t1 - t0, "\n"))
			t0 = t1
		annoEnt, idx = self.varanno(var)

		if debugOpt:
			t1 = [time.time()]
			sys.stderr.write("{}{}{}".format("BedAnno->anno [varanno] ... ", t1 - t0, "\n"))
			t0 = t1
		return annoEnt


	@staticmethod
	def _calfsOfst(trdbEnt, p):    # 0 based
		if "csta" not in trdbEnt or "nfs" not in trdbEnt:
			return 0

		ofst = 0
		for fsld in sorted(trdbEnt['nfs'].keys()):
			fsEd = int(fsld) + int(trdbEnt['nfs'][fsld])
			if int(fsld) <= int(p) and int(fsEd) < int(p):    # before p
				ofst =int(ofst)+ int(trdbEnt['nfs'][fsld])
			elif int(fsld) <= int(p) and int(fsEd) >= int(p):
				ofst =int(ofst) +int(fsld) + int(trdbEnt['nfs'][fsld]) - int(p)
				return [ofst]
			elif int(fsld) > int(p) and int(fsEd) <= int(p):
				return [ofst, (int(ofst) + int(trdbEnt['nfs'][fsld]))]
			else:
				break
		return [ofst]

	# check if region in duplicated frames are connected
	@staticmethod
	def _chkCodonPosRel(rsta, rsto):
		if rsta[2][0] == rsto[1][0] and rsta[2][1] <= rsto[1][1] + 1:
			# same codon connected frame
			return 1
		elif rsta[2][0] == rsto[1][0] + 1 and rsta[2][1] == 0 and rsto[1][1] == 2:
			# neighbor codon consecutive frame
			return 1
		elif rsta[2][0] < rsto[1][0]:
			# overlapped cases
			return 1
		else:
			return 0


	@staticmethod
	def _intronSO(r):
		if re.search(r'^I5U', r):
			return 'SO:0000447'
		elif re.search(r'^I3U', r):
			return 'SO:0000448'
		elif re.search(r'^I', r):
			return 'SO:0000191'
		else:
			raise "Not intron region tag [{}]".format(r)