#!/usr/bin/env python
# -*- coding:utf-8 -*-

from __future__ import division
import pysam, re,sys,csv,copy,random,logging,multiprocessing,os,functools
sys.dont_write_bytecode = True
#HG="/mnt/NL200/zhanglc/bin/BNC/program/NoahCare/db/alignment/hg19/hg19.fa"
#HG="/mnt/NL200/prod/repos/BNC/BNC/program/NoahCare/db/alignment/tgp_phase2_flat/hs37d5.fa"
debugging=False
f_ref,f_mut,f_N_both,f_N_consistent="","","",""
def LogExc(f):
    @functools.wraps(f)
    def wrapper(*largs, **kwargs):
        try:
            res=f(*largs, **kwargs)
        except Exception, e:
            multiprocessing.get_logger().exception(e)
            raise
        return res
    return wrapper
@LogExc
def mstch_count(a,b,hint):
	logger=multiprocessing.get_logger()
	mismatch_count = 0
	a_len = len(a)
	b_len = len(b)
	if (a_len != b_len):
		logger.error(hint)
		logger.error('error in mismatch_count: len(a)!=len(b)')
		logger.error("{}{}".format('a=', a))
		logger.error("{}{}".format('b=', b))
		sys.exit(1)
	for i in range(a_len):
		if(a[i] !=b [i]):
			mismatch_count += 1
	return mismatch_count


def get_ref(chr,start,end,HG):
	fastafile = pysam.Fastafile(HG)
	result=fastafile.fetch(chr, start, end)
	fastafile.close()
	return result

def is_dup(read):
	if(read.is_duplicate):
		return True
	else:
		return False


def get_le_de_ri_seq(chr,shift_start,shift_end,start,end,HG):
	flank_len=400
	l_seq=get_ref(chr, shift_start-flank_len, shift_start,HG)
	r_seq=get_ref(chr, shift_end, shift_end+flank_len,HG)
	del_ref_l=get_ref(chr, shift_end-flank_len, shift_end,HG)
	del_ref_r=get_ref(chr, shift_start, shift_start+flank_len,HG)
	l_del_ref=get_ref(chr,start-(flank_len-(shift_end-end)),start,HG)+get_ref(chr,end,shift_end,HG)
	r_del_ref=get_ref(chr,shift_start,start,HG)+get_ref(chr,end,end+(flank_len-(start-shift_start)),HG)
	return (l_seq,del_ref_l,del_ref_r,r_seq,l_del_ref,r_del_ref)


def find_only_one_position(read,seq):
	mismatch_dict=dict()
	for i in range(100):
		mismatch_dict[i]=list()
	for i in range(len(read)-len(seq)+1):
		mis_count=mstch_count(read[i:i+len(seq)],seq)
		mismatch_dict[mis_count].append(i)
	for i in range(100):
		if(len(mismatch_dict[i])==1):
			return True
		elif(len(mismatch_dict[i])>1):
			return False

def cigar_sum_index(cigar,read,del_parts,start,end,blocks,del_blocks,chr,HG):
	#print 'blocks,blocks,cigar,start,end=',blocks,cigar,start,end
	read_parts=list()
	del_part=list()
	D_index=list()
	I_index=list()
	M_index=list()
	S_index=list()
	cigar_pat=re.compile(r'[0-9]+D|[0-9]+M|[0-9]+S|[0-9]+I')
	D_M_array=cigar_pat.findall(cigar)
	DM_len=list()
	for item in D_M_array:
		DM_len.append(int(item[:-1]))  #['20M', '15D', '20M', '15D', '36M']
	sum=copy.deepcopy(DM_len)
	sum[0]=DM_len[0]
	for i in range(1,len(DM_len)):
		if("D" in D_M_array[i]):
			sum[i]=sum[i-1]
		else:
			sum[i]=sum[i-1]+DM_len[i]
	for i in range(len(D_M_array)):
		if "D" in D_M_array[i]:
			D_index.append(i)
		elif("I" in D_M_array[i]):
			I_index.append(i)
		elif("M" in D_M_array[i]):
			M_index.append(i)
		elif("S" in D_M_array[i]):
			S_index.append(i)
	cursor=0
	for i in range(len(D_M_array)):
		if("I" in D_M_array[i]):
			read_parts.append("")
		else:
			read_parts.append(read[cursor:sum[i]])
		cursor=sum[i]
	for i in range(len(D_M_array)):
		if("D" in D_M_array[i]):
			if(len(del_parts))>0:
				d_b_start,d_b_end=del_blocks[0]
				min_pos=min(start,end,d_b_start,d_b_end)
				max_pos=max(start,end,d_b_start,d_b_end)
				interval_min_max_ref=get_ref(chr,min_pos,max_pos,HG)
				ATGC_set=set()
				for j in interval_min_max_ref:
					ATGC_set.add(j)
				if(len(ATGC_set)>1):
					read_parts[i]=del_parts[0]
			else:
				break
			del_parts=del_parts[1:]
			del_blocks=del_blocks[1:]
	if("S" in cigar):
		for i in blocks:
			b_start,b_end=i
			if(b_start<=start and end<=b_end):
				for k in S_index:
					read_parts[k]=""
	new_read_start,new_read_end=blocks[0][0],blocks[-1][1]
	#if not (re.search(r'(\d+M\d+S)|(\d+S\d+M)',cigar) and (blocks[0][0] == start or blocks[-1][1]==end)):
	if(M_index[0]>0):
		for i in range(M_index[0])[::-1]:
			if("S" in D_M_array[i]):
				if(len(read_parts[i])>10) and new_read_start!=end: #####
					S_len=DM_len[i]
					S_add_seq=get_ref(chr,new_read_start-S_len,new_read_start,HG)
					read_parts[i]=S_add_seq
				new_read_start-=len(read_parts[i])
	if(M_index[-1]!=len(D_M_array)-1):
		for i in range(M_index[-1]+1,len(D_M_array)):
			if("S" in D_M_array[i]):
				if(len(read_parts[i])>10) and new_read_end!=start: ####
					S_len=DM_len[i]
					S_add_seq=get_ref(chr,new_read_end,new_read_end+S_len,HG)
					read_parts[i]=S_add_seq
				new_read_end+=len(read_parts[i])
	new_read="".join(read_parts)
	return new_read,new_read_start,new_read_end
	#else:
	#	return read,new_read_start,new_read_end

def insrt_DI_base_in_other_position(read,cigar,blocks,chr,start,end,shift_start,shift_end,HG):
	un_fill_in_D_len=0
	pre_start,pre_end=blocks[0]
	del_parts=list()
	del_blocks=list()
	for i in range(1,len(blocks)):
		cur_start,cur_end=blocks[i]
		if(shift_start<=pre_end and cur_start<=shift_end):
			del_blocks.append((pre_end,cur_start))
			del_parts.append("")
			un_fill_in_D_len=cur_start-pre_end
		else:
			if(pre_end!=cur_start):
				del_blocks.append((pre_end,cur_start))
				del_seq=get_ref(chr,pre_end,cur_start,HG)
				del_parts.append(del_seq)
		pre_start,pre_end=cur_start,cur_end
	new_read,new_read_start,new_read_end=cigar_sum_index(cigar,read,del_parts,start,end,blocks,del_blocks,chr,HG)
	return new_read,new_read_start,new_read_end,un_fill_in_D_len

def handle_repeat_dels(chr,start,end,HG):
	l_r_extend_100=get_ref(chr,start-100,end+100,HG)
	min_pos,max_pos=start,end
	l_start,l_end=start,end
	r_start,r_end=start,end
	d_len=end-start
	del_seq=get_ref(chr,start,end,HG)
	left_loop_mark=True
	right_loop_mark=True
	while(left_loop_mark or right_loop_mark):
		if(left_loop_mark):
			del_seq_l=get_ref(chr,l_start-d_len,l_end-d_len,HG)
			for i in range(d_len-1,-1,-1):
				if(del_seq[i]==del_seq_l[i]):
					min_pos-=1
				else:
					left_loop_mark=False
					break
				if(i==0):
					l_start-=d_len
					l_end-=d_len
		if(right_loop_mark):
			del_seq_r=get_ref(chr,r_start+d_len,r_end+d_len,HG)
			for i in range(d_len):
				if(del_seq[i]==del_seq_r[i]):
					max_pos+=1
				else:
					right_loop_mark=False
					break
				if(i==(d_len-1)):
					r_start+=d_len
					r_end+=d_len
	return (min_pos,max_pos)

def is_interset_le_shift_interval(shift_start,shift_end,read_start,read_end,read_seq,HG,chr,c_gar,start,end):
	max_value=max(start,read_start)
	min_value=min(end,read_end)
	if(max_value<min_value):
		if re.search(r'^\d+M$',c_gar):
			shift_inv_seq=""
			ref_seq=""
			if max_value==start and read_end>start and read_end<=end:
				shift_inv_seq=read_seq[-(read_end-start):]
				m=re.search(r'(N+)$',shift_inv_seq)
				if m:
					shift_inv_seq=shift_inv_seq[:-len(m.group(1))]
					ref_seq=get_ref(chr,max_value,min_value-len(m.group(1)),HG)
				else:
					ref_seq=get_ref(chr,max_value,min_value,HG)
			elif max_value==read_start and read_start>=start and read_start<end:
				shift_inv_seq=read_seq[:(end-read_start)]
				m=re.search(r'^(N+)',shift_inv_seq)
				if m:
					shift_inv_seq=shift_inv_seq[len(m.group(1)):]
					ref_seq=get_ref(chr,max_value+len(m.group(1)),min_value,HG)
				else:
					ref_seq=get_ref(chr,max_value,min_value,HG)
			if shift_inv_seq!=ref_seq:
				return False
	
	shift_len=shift_end-shift_start
	max_value=max(shift_start,read_start)
	min_value=min(shift_end,read_end)
	if(max_value<min_value):
		'''
		if re.search(r'^\d+M$',c_gar):
			shift_inv_seq=""
			ref_seq=""
			if max_value==shift_start:
				shift_inv_seq=read_seq[:shift_end-read_start]
				m=re.search(r'^(N+)')
				if m:
					shift_inv_seq=shift_inv_seq[len(m.group(1)):]
					ref_seq=get_ref(chr,max_value+len(m.group(1)),min_value,HG)
				else:
					ref_seq=get_ref(chr,max_value,min_value,HG)
			elif max_value==read_start:
				shift_inv_seq=read_seq[-(read_end-shift_start):]
				m=re.search(r'(N+)$')
				if m:
					shift_inv_seq=shift_inv_seq[:-len(m.group(1))]
					ref_seq=get_ref(chr,max_value,min_value-len(m.group(1)),HG)
				else:
					ref_seq=get_ref(chr,max_value,min_value,HG)
			if shift_inv_seq!=ref_seq:
				return False
		'''
		intersect_len=min_value-max_value
		if(intersect_len<shift_len):
			return True
	else:
		return False



def get_seqs(start,end,shift_start,shift_end,M_start,M_end,new_read,un_fill_in_D_len):
	left_len=shift_start-M_start
	right_len=M_end-shift_end
	l_before_shift_end=new_read[:(shift_end-M_start-un_fill_in_D_len)]
	r_after_shift_start=new_read[-(M_end-shift_start-un_fill_in_D_len):]
	return l_before_shift_end,r_after_shift_start
@LogExc
def is_mut(normal_ref_left,normal_ref_right,del_ref_left,del_ref_right,l_before_shift_end,r_after_shift_start,c_gar,blocks,read,chr,start,end):
	flag_other=0
	l_r_mismatch_count=None
	
	l_r_mismatch_rate=None
	def print_variables(f_handle):
		savedStdout = sys.stdout
		sys.stdout=f_handle
		print("read=",read)
		print("l_before_shift_end=",l_before_shift_end)
		print("normal_ref_left=   ",normal_ref_left)
		print("del_ref_left=      ",del_ref_left)
		print("left_normal_mismatch=",left_normal_mismatch)
		print("left_del_mismatch=",left_del_mismatch)
		print("r_after_shift_start=",r_after_shift_start)
		print("normal_ref_right=   ",normal_ref_right)
		print("del_ref_right=      ",del_ref_right)
		print("right_normal_mismatch=",right_normal_mismatch)
		print("right_del_mismatch=",right_del_mismatch)
		sys.stdout = savedStdout
	normal_ref_left=normal_ref_left[-len(l_before_shift_end):]
	normal_ref_right=normal_ref_right[:len(r_after_shift_start)]
	del_ref_left=del_ref_left[-len(l_before_shift_end):]
	del_ref_right=del_ref_right[:len(r_after_shift_start)]
	left_del_mismatch=mstch_count(l_before_shift_end,del_ref_left,"304")
	left_normal_mismatch=mstch_count(l_before_shift_end,normal_ref_left,"305")
	right_del_mismatch=mstch_count(r_after_shift_start,del_ref_right,"306")
	right_normal_mismatch=mstch_count(r_after_shift_start,normal_ref_right,"307")
	if re.search(r'^\d+M$',c_gar) and blocks[0][1]>start and blocks[0][1]<=end:
		left_del_mismatch=0
		left_normal_mismatch=0
	if re.search(r'^\d+M$',c_gar) and blocks[0][0]<end and blocks[0][0]>=start:
		right_del_mismatch=0
		right_normal_mismatch=0
	if re.search(r'^\d+S\d+M$',c_gar) and blocks[0][0]==end:
		right_del_mismatch=0
		right_normal_mismatch=0
	if re.search(r'^\d+M\d+S$',c_gar) and blocks[-1][1]==start:
		left_del_mismatch=0
		left_normal_mismatch=0
	l_r_read_len=len(l_before_shift_end)+len(r_after_shift_start)
	print_variables(sys.stdout) if debugging else False
	f_out_handle,return_result="",""
	if(left_del_mismatch<left_normal_mismatch and right_del_mismatch<right_normal_mismatch): #and right_del_mismatch<right_normal_mismatch
		l_r_mismatch_count=left_del_mismatch+right_del_mismatch
		f_out_handle=f_mut
		return_result=1
	elif(left_del_mismatch>left_normal_mismatch and right_del_mismatch>right_normal_mismatch): #and right_del_mismatch>right_normal_mismatch
		l_r_mismatch_count=left_normal_mismatch+right_normal_mismatch
		f_out_handle=f_ref
		return_result=0
	elif(left_del_mismatch>left_normal_mismatch and right_del_mismatch<right_normal_mismatch):
		if re.search(r"^\d+M$", c_gar):

			l_r_mismatch_count=0
			f_out_handle=f_ref
			return_result=0
		else:
			if left_del_mismatch+right_del_mismatch<left_normal_mismatch+right_normal_mismatch:
				l_r_mismatch_count=left_del_mismatch+right_del_mismatch
				f_out_handle=f_mut
				return_result=1
			elif left_del_mismatch+right_del_mismatch>left_normal_mismatch+right_normal_mismatch:
				l_r_mismatch_count=left_normal_mismatch+right_normal_mismatch
				f_out_handle=f_ref
				return_result=0
			else:
				f_out_handle=f_N_consistent
				return_result=-1
				flag_other=1
	elif(left_del_mismatch<left_normal_mismatch and right_del_mismatch>right_normal_mismatch):
		#if("M" in c_gar and "D" not in c_gar):
		if re.search(r"^\d+M$", c_gar):

			l_r_mismatch_count=0
			f_out_handle=f_ref
			return_result=0
		else:
			if left_del_mismatch+right_del_mismatch<left_normal_mismatch+right_normal_mismatch:
				l_r_mismatch_count=left_del_mismatch+right_del_mismatch
				f_out_handle=f_mut
				return_result=1
			elif left_del_mismatch+right_del_mismatch>left_normal_mismatch+right_normal_mismatch:
				l_r_mismatch_count=left_normal_mismatch+right_normal_mismatch
				f_out_handle=f_ref
				return_result=0
			else:
				f_out_handle=f_N_consistent
				return_result=-1
				flag_other=1
	elif(left_del_mismatch==left_normal_mismatch and right_del_mismatch<right_normal_mismatch):
		l_r_mismatch_count=left_del_mismatch+right_del_mismatch
		f_out_handle=f_mut
		return_result=1
	elif(left_del_mismatch==left_normal_mismatch and right_del_mismatch>right_normal_mismatch):
		l_r_mismatch_count=left_del_mismatch+right_normal_mismatch
		f_out_handle=f_ref
		return_result=0
	elif(left_del_mismatch>left_normal_mismatch and right_del_mismatch==right_normal_mismatch):
		l_r_mismatch_count=left_normal_mismatch+right_del_mismatch
		f_out_handle=f_ref
		return_result=0
	elif(left_del_mismatch<left_normal_mismatch and right_del_mismatch==right_normal_mismatch):
		l_r_mismatch_count=left_del_mismatch+right_del_mismatch
		f_out_handle=f_mut
		return_result=1
	else:
		f_out_handle=f_N_both
		return_result=-2
		flag_other=1
	print_variables(f_out_handle) if debugging else False
	if return_result==0 or return_result==1:
		l_r_mismatch_rate=(l_r_mismatch_count+0.0)/l_r_read_len
	return return_result,l_r_mismatch_rate,flag_other


def is_not_exact_D_len(cigar,start,end,blocks):
	m=re.search(r'^([0-9]+)M([0-9]+)D([0-9]+)M$',cigar)
	if m and int(blocks[0][1])==start and int(blocks[1][0])==end and int(m.group(2))==(end-start):
		return False
	else:
		return True
	'''
	D_len=end-start
	cigar_pat=re.compile(r'[0-9]+D')
	D_array=cigar_pat.findall(cigar)
	for i in D_array:
		lens=int(i[:-1])
		if(lens==D_len):
			return False
	return True
	'''


@LogExc
def correct_under_del_rate(bam,chr,start,end,ref_genome):
	ref_qname_mismatch_dict=dict()
	alt_qname_mismatch_dict=dict()
	ref_query_name_set=set()
	alt_query_name_set=set()
	HG=ref_genome
	all_reads=0
	samfile = pysam.AlignmentFile(bam, "rb")
	ref_count=0
	mut_count=0
	N_mut_and_ref=0
	inconsistent_count=0
	line_num=0
	d_seq=get_ref(chr,start,end,HG)
	left_before_d_seq=get_ref(chr,start-100,start,HG)
	right_after_d_seq=get_ref(chr,end,end+100,HG)
	shift_start,shift_end=handle_repeat_dels(chr,start,end,HG)
	del_seq_left_once,del_seq_l_once,del_seq_r_once,del_seq_right_once,l_del_ref,r_del_ref=get_le_de_ri_seq(chr, shift_start, shift_end,start,end,HG)
	del_seq_left_once_compare,del_seq_right_once_compare="",""
	soft_index=0
	other_reads_count=0
	other_reads_count_set=set()
	other_reads_count_qname_set=set()
	'''
	reads_other="{}-{}-{}-{}-other.bam".format(chr,start,end,d_seq)
	reads_ref="{}-{}-{}-{}-ref.bam".format(chr,start,end,d_seq)
	reads_alt="{}-{}-{}-{}-alt.bam".format(chr,start,end,d_seq)
	reads_other_f = pysam.AlignmentFile(reads_other, "wb", template=samfile)
	reads_ref_f = pysam.AlignmentFile(reads_ref, "wb", template=samfile)
	reads_alt_f = pysam.AlignmentFile(reads_alt, "wb", template=samfile)
	'''
	for read in samfile.fetch(chr, start-30, end+30):
		def print_variables(f_handle):
			savedStdout = sys.stdout
			sys.stdout=f_handle
			print('c_gar=',c_gar)
			print('d_seq=',d_seq)
			print('query_name=',query_name)
			print("read_seq=",read_seq)
			print("blocks=",blocks)
			print("shift_start=",shift_start)
			print("shift_end=",shift_end)
			print("read_sequence_ori=",read_sequence_ori)
			print("read_sequence_ori_start=",read_sequence_ori_start)
			print("read_sequence_ori_end=",read_sequence_ori_end)
			print("ref_count=",ref_count)
			print("mut_count=",mut_count)
			print("N_mut_and_ref=",N_mut_and_ref)
			print("*"*50)
			sys.stdout = savedStdout
		query_name=read.query_name
		c_gar=read.cigarstring
		if(read.is_unmapped or read.mate_is_unmapped):
			continue
		blocks=read.get_blocks() 
		if(is_dup(read)):
			continue
		if "M" not in c_gar:
			continue
		if read.reference_id!=read.next_reference_id:
			continue
		if "D" in c_gar:
			indictator=0
			flag=False
			for i in range(len(blocks)-1):
				now_start,now_end=blocks[i][0],blocks[i][1]
				next_start,next_end=blocks[i+1][0],blocks[i+1][1]
				if now_end==start and next_start!=end and now_end!=next_start:
					indictator=1
					break
				if now_end==start and next_start==end and now_end!=next_start:
					flag=True
					break
			if flag:
				alt_query_name_set.add(query_name)
				mut_count+=1
				alt_qname_mismatch_dict[query_name]=-1
				#reads_alt_f.write(read)
				continue
			if indictator:
				#reads_other_f.write(read)
				other_reads_count_qname_set.add(query_name)
				other_reads_count_set.add(query_name)
				other_reads_count+=1
				continue
		
		cigar_s=read.cigarstring 
		read_seq=read.query_sequence 
		read_sequence_ori=read_seq
		'''
		read_sequence_ori_start,read_sequence_ori_end=0,0
		read_sequence_ori,read_sequence_ori_start,read_sequence_ori_end,un_fill_in_D_len=insrt_DI_base_in_other_position(read_sequence_ori,c_gar,blocks,chr,start,end,shift_start,shift_end,HG)
		'''
		if re.search(r'^\d+M$',c_gar) and (blocks[0][0]>=end or blocks[0][1]<=start ):
			continue
		'''
		if re.search(r'^\d+M$',c_gar) and "N" in read_seq:
			r_N_nb, l_N_nb="", ""
			m=re.search(r'(N+)$',read_seq)
			if m:
				r_N_nb=m.group(1)
				if blocks[0][1]-len(r_N_nb)<=start:
					continue
			m=re.search(r'^(N+)',read_seq)
			if m:
				l_N_nb=m.group(1)
				if blocks[0][0]+len(l_N_nb)>=end:
					continue
		'''
		m1=re.search(r'(\d+)M$',c_gar)
		m2=re.search(r'^(\d+)M',c_gar)
		if (m1 or m2) and "N" in read_seq:
			if m1:
				M_seq=read_seq[-int(m1.group(1)):]
				m=re.search(r'(N+)$',M_seq)
				if m:
					r_N_nb=m.group(1)
					if blocks[-1][1]-len(r_N_nb)<=start:
						continue
			if m2:
				M_seq=read_seq[:int(m2.group(1))]
				m=re.search(r'^(N+)',read_seq)
				if m:
					l_N_nb=m.group(1)
					if blocks[0][0]+len(l_N_nb)>=end:
						continue
			#r_N_nb, l_N_nb="", ""
		
		'''
		if re.search(r'^\d+S\d+M\d+S$',c_gar):
			min_max=max(blocks[0][1],end)
			max_min=min(start,blocks[0][0])
			if min_max>max_min:
				if blocks[0][0]!=end or blocks[0][1]!=start:
					reads_ref_f.write(read)
					ref_query_name_set.add(query_name)
					ref_qname_mismatch_dict[query_name]=-1
					continue
		''' 
		if re.search(r'^\d+M\d+H$|^\d+H\d+M$',c_gar):
			min_max=max(blocks[0][1],end)
			max_min=min(start,blocks[0][0])
			if min_max>max_min:
				#reads_ref_f.write(read)
				ref_query_name_set.add(query_name)
				ref_qname_mismatch_dict[query_name]=-1
				continue
			else:
				continue
		all_reads+=1
		#if re.search(r'(^\d+M\d+S$)|(^\d+S\d+M$)',c_gar) :
		
		m= re.search(r'^(\d+)S(\d+)M\d*',c_gar)
		if m:
			n=re.search(r'^(N+)',read_seq)
			N_nr=0
			if n:
				N_nr=len(n.group(1))
			S_len=int(m.group(1))-N_nr
			if blocks[0][0]==end:
				if S_len==0:
					continue
				if read_seq[N_nr:int(m.group(1))]==left_before_d_seq[-S_len:]:
					soft_index+=1
					alt_query_name_set.add(query_name)
					mut_count+=1
					alt_qname_mismatch_dict[query_name]=-1
					#reads_alt_f.write(read)
					continue
				msmch=mstch_count(read_seq[N_nr:int(m.group(1))],left_before_d_seq[-S_len:],"560")
				if (S_len>=3 and msmch>=3) or  msmch/S_len>=0.5:
					#reads_ref_f.write(read)
					ref_query_name_set.add(query_name)
					ref_qname_mismatch_dict[query_name]=-1
					continue
			else:
				if blocks[0][0]>end and blocks[0][0]-S_len>=end:
					continue
				if re.search(r'(\d+)M$',c_gar) and blocks[-1][1]<=start:
					continue
				#if (blocks[0][0]>end and blocks[0][0]-S_len<end) or blocks[0][0]<end: 
				if (blocks[0][0]>end and blocks[0][0]-S_len<end) : 
					#reads_ref_f.write(read)
					ref_query_name_set.add(query_name)
					ref_qname_mismatch_dict[query_name]=-1
					continue
				if blocks[0][0]<end:
					r_most_pos=blocks[-1][1]
					q=re.search(r'(\d+)S$',c_gar)
					if q:
						n=re.search(r'(N+)$',read_seq)
						N_nr=0
						if n:
							N_nr=len(n.group(1))
						r_most_pos+=int(q.group(1))-N_nr
					if r_most_pos>start:
						#reads_ref_f.write(read)
						ref_query_name_set.add(query_name)
						ref_qname_mismatch_dict[query_name]=-1
						continue
					else:
						continue
				if S_len==0 and blocks[0][0]>=end:
					continue
		
		m= re.search(r'\d*[A-Z]*(\d+)M(\d+)S$',c_gar)
		if m :
			n=re.search(r'(N+)$',read_seq)
			N_nr=0
			if n:
				N_nr=len(n.group(1))
			S_len=int(m.group(2))-N_nr
			if blocks[-1][1]==start:
				if S_len==0:
					continue
				if read_seq[-int(m.group(2)):-N_nr]==right_after_d_seq[:int(m.group(2))-N_nr]:
					soft_index+=1
					alt_query_name_set.add(query_name)
					alt_qname_mismatch_dict[query_name]=-1
					#reads_alt_f.write(read)
					continue
				msmch=0
				if N_nr==0:
					msmch=mstch_count(read_seq[-int(m.group(2)):],right_after_d_seq[:int(m.group(2))],"599")
				else:
					msmch=mstch_count(read_seq[-int(m.group(2)):-N_nr],right_after_d_seq[:int(m.group(2))-N_nr],"601")
				if (S_len>=3 and msmch>=3) or  msmch/S_len>=0.5:
					#reads_ref_f.write(read)
					ref_query_name_set.add(query_name)
					ref_qname_mismatch_dict[query_name]=-1
					continue
			else:
				if blocks[-1][1]<start and blocks[-1][1]+S_len<=start:
					continue
				if re.search(r'^(\d+)M',c_gar) and blocks[0][0]>=end:
					continue
				#if (blocks[-1][1] < start and blocks[-1][1]+S_len>start) or blocks[-1][1]>start:
				if (blocks[-1][1] < start and blocks[-1][1]+S_len>start) :
					#reads_ref_f.write(read)
					ref_query_name_set.add(query_name)
					ref_qname_mismatch_dict[query_name]=-1
					continue
				if blocks[-1][1]>start:
					l_most_pos=blocks[0][0]
					q=re.search(r'^(\d+)S',c_gar)
					if q:
						n=re.search(r'^(N+)',read_seq)
						N_nr=0
						if n:
							N_nr=len(n.group(1))
						l_most_pos-=int(q.group(1))-N_nr
					if l_most_pos<end:
						#reads_ref_f.write(read)
						ref_query_name_set.add(query_name)
						ref_qname_mismatch_dict[query_name]=-1
						continue
					else:
						continue
				if S_len==0 and blocks[-1][1]<=start:
					continue
				
		'''
		m=re.search(r'^(\d+)S\d+M\d+',c_gar)
		if m :
			n=re.search(r'^(N+)',read_seq)
			N_nr=0
			if n:
				N_nr=len(n.group(1))
			S_len=int(m.group(1))-N_nr
			if blocks[0][0]==end:
				if S_len==0:
					continue
				if read_seq[N_nr:int(m.group(1))]==left_before_d_seq[-S_len:]:
					alt_query_name_set.add(query_name)
					mut_count+=1
					alt_qname_mismatch_dict[query_name]=-1
					reads_alt_f.write(read)
					continue
				msmch=mstch_count(read_seq[N_nr:int(m.group(1))],left_before_d_seq[-S_len:],"631")
				if (int(m.group(1))>=3 and msmch>=3) or  msmch/int(m.group(1))>=0.5:
					reads_ref_f.write(read)
					ref_query_name_set.add(query_name)
					ref_qname_mismatch_dict[query_name]=-1
					continue
			else:
				if blocks[0][0]>end and blocks[0][0]-S_len<end:
					reads_ref_f.write(read)
					ref_query_name_set.add(query_name)
					ref_qname_mismatch_dict[query_name]=-1
					continue
				if S_len==0 and blocks[0][0]>=end:
					continue
		m=re.search(r'\d+[A-Z]\d+M(\d+)S$',c_gar)
		if m :
			n=re.search(r'(N+)$',read_seq)
			N_nr=0
			if n:
				N_nr=len(n.group(1))
			S_len=int(m.group(1))-N_nr
			if blocks[-1][1]==start:
				if S_len==0:
					continue
				if read_seq[-int(m.group(1)):-N_nr]==right_after_d_seq[:S_len]:
					alt_query_name_set.add(query_name)
					mut_count+=1
					alt_qname_mismatch_dict[query_name]=-1
					reads_alt_f.write(read)
					continue
				msmch=0
				if N_nr==0:
					msmch=mstch_count(read_seq[-int(m.group(1)):],right_after_d_seq[:S_len],"663")
				else:
					msmch=mstch_count(read_seq[-int(m.group(1)):-N_nr],right_after_d_seq[:S_len],"665")
				if (int(m.group(1))>=3 and msmch>=3) or  msmch/int(m.group(1))>=0.5:
					reads_ref_f.write(read)
					ref_query_name_set.add(query_name)
					ref_qname_mismatch_dict[query_name]=-1
					continue
			else:
				if blocks[-1][1]<start and blocks[-1][1]+S_len>start:
					reads_ref_f.write(read)
					ref_query_name_set.add(query_name)
					ref_qname_mismatch_dict[query_name]=-1
					continue
				if S_len==0 and blocks[-1][1]<=start:
					continue
		'''
		read_sequence_ori_start,read_sequence_ori_end=0,0
		read_sequence_ori,read_sequence_ori_start,read_sequence_ori_end,un_fill_in_D_len=insrt_DI_base_in_other_position(read_sequence_ori,c_gar,blocks,chr,start,end,shift_start,shift_end,HG)
		if read_sequence_ori_start >=end or read_sequence_ori_end<=start:
			continue
		
		if not(start==shift_start and end==shift_end):
			next_or_coni=is_interset_le_shift_interval(shift_start,shift_end,read_sequence_ori_start,read_sequence_ori_end,read_sequence_ori,HG,chr,c_gar,start,end)
			if(next_or_coni):
				#reads_other_f.write(read)
				other_reads_count_qname_set.add(query_name)
				other_reads_count+=1
				other_reads_count_set.add(query_name)
				continue
		line_num+=1
		read_sequence_ori_len=len(read_sequence_ori)
		l_before_shift_end,r_after_shift_start=get_seqs(start,end,shift_start,shift_end,read_sequence_ori_start,read_sequence_ori_end,read_sequence_ori,un_fill_in_D_len)
		m=re.search(r'^([0-9]+)M([0-9]+)D([0-9]+)M$',c_gar)
		if m :
			if(is_not_exact_D_len(c_gar,start,end,blocks)):
				#reads_other_f.write(read)
				other_reads_count_qname_set.add(query_name)
				other_reads_count+=1
				other_reads_count_set.add(query_name)
				continue
		m1=re.search(r'^\d+S(\d+)M(\d+)D(\d+)M$',c_gar)
		m2=re.search(r'^(\d+)M(\d+)D(\d+)M\d+S$',c_gar)
		MDM_D_len=0
		if m1:
			MDM_D_len=int(m1.group(2))
		if m2:
			MDM_D_len=int(m2.group(2))
		if m1 or m2:
			if blocks[0][1]==start or blocks[1][0]==end or (blocks[0][1]>start and blocks[1][0]<end):
				if MDM_D_len!=end-start:
					#reads_ref_f.write(read)
					ref_query_name_set.add(query_name)
					ref_qname_mismatch_dict[query_name]=-1
					continue
		mut_ref_N_bool,mismatch_rate=None,None
		mut_ref_N_bool,mismatch_rate,flag_other=is_mut(del_seq_l_once,del_seq_r_once,l_del_ref,r_del_ref,l_before_shift_end,r_after_shift_start,c_gar,blocks,read_sequence_ori,chr,start,end)
		other_reads_count+=flag_other
		if flag_other==1:
			#reads_other_f.write(read)
			other_reads_count_qname_set.add(query_name)
		f_out_handle=""
		if(mut_ref_N_bool==1):
			alt_qname_mismatch_dict[query_name]=mismatch_rate
			alt_query_name_set.add(query_name)
			mut_count+=1
			f_out_handle=f_mut
			#reads_alt_f.write(read)
		elif(mut_ref_N_bool==0):
			ref_qname_mismatch_dict[query_name]=mismatch_rate
			ref_query_name_set.add(query_name)
			ref_count+=1
			f_out_handle=f_ref
			#reads_ref_f.write(read)
		elif(mut_ref_N_bool==-1):
			inconsistent_count+=1
			f_out_handle=f_N_consistent
		elif(mut_ref_N_bool==-2):
			N_mut_and_ref+=1
			f_out_handle=f_N_both
		print_variables(f_out_handle) if debugging else False
		print_variables(sys.stdout) if debugging else False
	samfile.close() 
	'''
	reads_other_f.close()
	reads_ref_f.close()
	reads_alt_f.close()
	
	os.system('samtools index '+reads_other +" "+reads_other+".bai")
	os.system('samtools index '+reads_ref +" "+reads_ref+".bai")
	os.system('samtools index '+reads_alt +" "+reads_alt+".bai")
	'''
	for set_item in ref_query_name_set & alt_query_name_set:
		if ref_qname_mismatch_dict[set_item] < alt_qname_mismatch_dict[set_item]:
			alt_query_name_set.remove(set_item)
		elif ref_qname_mismatch_dict[set_item] > alt_qname_mismatch_dict[set_item]:
			ref_query_name_set.remove(set_item)
		else:
			random_data=random.randint(1,2)
			if random_data==1:
				alt_query_name_set.remove(set_item)
			else:
				ref_query_name_set.remove(set_item)
	
	if(len(alt_query_name_set)+len(ref_query_name_set)+other_reads_count==0):
		result=0
	else:
		result=len(alt_query_name_set)/(len(alt_query_name_set)+len(ref_query_name_set)+other_reads_count)
	#for ii in alt_query_name_set:
	#	print ii
	if(debugging):
		print('filter line_num=',line_num)
		print("filter all_reads=",all_reads)
		print("filter mut_count=",mut_count)
		print("filter ref_count=",ref_count)
		print("filter N_mut_and_ref=",N_mut_and_ref)
		print("filter mut_count/(mut_count+ref_count)=",result)
	'''
	file="alt_read_name.txt"
	Name_H=open(file,"w")
	for i in alt_query_name_set:
		print>>Name_H,i
	Name_H.close()
	'''
	return result,len(alt_query_name_set),len(alt_query_name_set)+len(ref_query_name_set)+other_reads_count



if __name__ == "__main__":
	log_file="test_debug"
	f=logging.Formatter('[%(levelname)s %(processName)s %(asctime)s %(funcName)s] %(message)s')
	h=logging.FileHandler(log_file, 'w')
	h.setFormatter(f)
	debugging=False
	csv_file="160005583BCD_160005583BPD.snv_idl.anno.csv"
	ref_genome="/mnt/NL200/prod/repos/BNC/BNC/program/NoahCare/db/alignment/tgp_phase2_flat/hs37d5.fa"
	#ref_genome="/mnt/NL200/refgenomes/tmp_hg19/hg19.fa"
	
	bam="/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000559FD/output/BACKUP/Bamdir/170419_BS1001160019_CL100015943_L001_V1D5HUMRNETKPPEI-411-v14.sort.markdup.realign.recald.bam"
	bam0="""
	/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000558FD/output/BACKUP/Bamdir/170419_BS1001160019_CL100015943_L001_V1D5HUMRNETKOPEI-410-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000559FD/output/BACKUP/Bamdir/170419_BS1001160019_CL100015943_L001_V1D5HUMRNETKPPEI-411-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000568FD/output/BACKUP/Bamdir/170419_BS1001160019_CL100015943_L001_V1D5HUMRNETKXPEI-420-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000571FD/output/BACKUP/Bamdir/170419_BS1001160019_CL100015943_L001_V1D5HUMRNETKZPEI-422-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000601FD/output/BACKUP/Bamdir/170419_BS1001160019_CL100015943_L001_V1D5HUMRNETMCPEI-406-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000603FD/output/BACKUP/Bamdir/170421_BS1001160015_CL100015950_L001_V1D5HUMRNETMDPEI-402-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000605FD/output/BACKUP/Bamdir/170421_BS1001160015_CL100015950_L001_V1D5HUMRNETQCPEI-404-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000606FD/output/BACKUP/Bamdir/170421_BS1001160015_CL100015950_L001_V1D5HUMRNETQDPEI-405-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000608FD/output/BACKUP/Bamdir/170421_BS1001160015_CL100015950_L001_V1D5HUMRNETQFPEI-407-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000627FD/output/BACKUP/Bamdir/170421_BS1001160015_CL100015950_L001_V1D5HUMRNETQVPEI-424-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000648FD/output/BACKUP/Bamdir/170421_BS1001160015_CL100015950_L001_V1D5HUMRNETRJPEI-445-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000649FD/output/BACKUP/Bamdir/170421_BS1001160015_CL100015950_L001_V1D5HUMRNETRKPEI-446-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000668FD/output/BACKUP/Bamdir/170629_BS1001160241_CL100021692_L001_V1D5HUMRNFAUDPEI-413-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000671FD/output/BACKUP/Bamdir/170629_BS1001160241_CL100021378_L001_V1D5HUMRNFAUWPEI-434-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000674FD/output/BACKUP/Bamdir/170629_BS1001160241_CL100021378_L001_V1D5HUMRNFAVPPEI-405-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000676FD/output/BACKUP/Bamdir/170629_BS1001160241_CL100021378_L001_V1D5HUMRNFAVFPEI-443-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000677FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVKQPEI-419-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000680FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVKRPEI-420-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000691FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVKLPEI-414-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000695FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVKYPEI-427-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000702FD/output/BACKUP/Bamdir/170513_BS1001160015_CL100015948_L001_V1D5HUMRNEVKVPEI-424-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000704FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVLAPEI-429-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000748FD/output/BACKUP/Bamdir/170513_BS1001160015_CL100016039_L001_V1D5HUMRNEWDAPEI-426-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001023FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVKKPEI-413-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001042FD/output/BACKUP/Bamdir/170513_BS1001160015_CL100016039_L001_V1D5HUMRNEWBIPEI-430-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001058FD/output/BACKUP/Bamdir/170513_BS1001160015_CL100015948_L001_V1D5HUMRNEVZNPEI-431-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001060FD/output/BACKUP/Bamdir/170513_BS1001160015_CL100016039_L001_V1D5HUMRNEVKAPEI-403-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001063FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVIOPEI-413-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001064FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVJJPEI-434-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001067FD/output/BACKUP/Bamdir/170513_BS1001160015_CL100016039_L001_V1D5HUMRNEWCYPEI-424-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001070FD/output/BACKUP/Bamdir/170513_BS1001160015_CL100015948_L001_V1D5HUMRNEWBLPEI-433-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001074FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVIEPEI-403-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001089FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVIDPEI-402-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001105FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVJIPEI-433-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001153FD/output/BACKUP/Bamdir/170513_BS1001160015_CL100016039_L001_V1D5HUMRNEWCCPEI-402-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001169FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVJGPEI-431-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001215FD/output/BACKUP/Bamdir/170629_BS1001160241_CL100021692_L001_V1D5HUMRNFATDPEI-435-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160009665FD/output/BACKUP/Bamdir/170221_NS500823_H53WJBGX2_L001_V1D5HUMCNEOMXPEI-211-v.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170005173FD/output/BACKUP/Bamdir/170422_TPNB500129_HKKKFBGX2_L001_HUMCNEULVPEI-263-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160007033FD/output/BACKUP/Bamdir/161104_TPNB500129_HV3T2BGXY_L001_V1D5HUMCAEFUYPEI-253-v.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160007559FD/output/BACKUP/Bamdir/161108_NS500823_H33L5BGX2_L001_V1D5HUMCAEGCAPEI-213-v.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160007942FD/output/BACKUP/Bamdir/161202_NS500823_H7NW7BGX2_L001_V1D5HUMCAEIXGPEI-252-v.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160007603FD/output/BACKUP/Bamdir/161212_NS500823_H7KNVBGX2_L001_V1D5HUMCAEJUUPEI-202-v13.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160009409FD/output/BACKUP/Bamdir/161220_TPNB500129_H7NHHBGX2_L001_V1D5HUMCAEKEEPEI-241-v106.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160004217FD/output/BACKUP/Bamdir/161224_TPNB500129_HF27LBGX2_L001_V1D5HUMCAEKNGPEI-239-v.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160007421FD/output/BACKUP/Bamdir/161231_TPNB500129_HF2HFBGX2_L001_V1D5HUMCAEKXNPEI-203-v106.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160009690FD/output/BACKUP/Bamdir/170103_TPNB500129_HCYW3BGX2_L001_V1D5HUMCAEKZXPEI-216-v106.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160007236FD/output/BACKUP/Bamdir/170109_TPNB500129_HCYNKBGX2_L001_V1D5HUMCAELNTPEI-259-v106.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160002650FD/output/BACKUP/Bamdir/170110_NS500823_HF2VWBGX2_L001_V1D5HUMCAELOLPEI-269-v.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160005892FD/output/BACKUP/Bamdir/170208_NS500823_H7YW3AFXX_L001_V1D5HUMCNENGQPEI-203-v106.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160009010FD/output/BACKUP/Bamdir/170211_NS500823_H5MLCBGX2_L001_V1D5HUMCNENPUPEI-219-v.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160002663FD/output/BACKUP/Bamdir/170218_NS500823_H53L2BGX2_L001_V1D5HUMCNEOJLPEI-224-v106.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170000362FD/output/BACKUP/Bamdir/170224_NS500823_H57CNBGX2_L001_V1D5HUMCNEOVWPEI-253-v.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160009403FD/output/BACKUP/Bamdir/170225_NS500823_H5MFYBGX2_L001_V1D5HUMCNEOYSPEI-295-v106.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160009149FD/output/BACKUP/Bamdir/170226_TPNB500129_H53F5BGX2_L001_V1D5HUMCNEOANPEI-228-v.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160005885FD/output/BACKUP/Bamdir/170312_TPNB500129_H5K5FBGX2_L001_HUMCNEQPPPEI-247-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160002638FD/output/BACKUP/Bamdir/170321_TPNB500129_H5NHKBGX2_L001_HUMCNERLCPEI-227-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170002320FD/output/BACKUP/Bamdir/170325_TPNB500129_H5MTTBGX2_L001_V1D5HUMCNERWWPEI-261-v13.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170005652FD/output/BACKUP/Bamdir/170330_TPNB500129_H5MJJBGX2_L001_V1D5HUMCNESJDPEI-285-v13.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170005864FD/output/BACKUP/Bamdir/170403_TPNB500129_H5CKCBGX2_L001_HUMCNESUZPEI-221-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170005060FD/output/BACKUP/Bamdir/170408_TPNB500129_HKCHFBGX2_L001_HUMCNETBVPEI-287-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170006964FD/output/BACKUP/Bamdir/170410_TPNB500129_HKHCMBGX2_L001_HUMCNETEIPEI-202-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170005061FD/output/BACKUP/Bamdir/170411_TPNB500129_HKHJMBGX2_L001_V1D5HUMCNETGVPEI-220-v13.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170002604FD/output/BACKUP/Bamdir/170413_TPNB500129_HKHKJBGX2_L001_HUMCNETNEPEI-272-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170005063FD/output/BACKUP/Bamdir/170413_TPNB500129_HKHKJBGX2_L001_V1D5HUMCNETNDPEI-255-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170005400FD/output/BACKUP/Bamdir/170421_TPNB500129_HKHJKBGX2_L001_HUMCNEUIHPEI-293-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170006773FD/output/BACKUP/Bamdir/170421_TPNB500129_HKHJKBGX2_L001_HUMCNEUIGPEI-218-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170005156FD/output/BACKUP/Bamdir/170426_TPNB500129_HJV7HBGX2_L001_HUMCNEUVUPEI-295-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170006909FD/output/BACKUP/Bamdir/170428_TPNB500129_HVCNJBGX2_L001_HUMCNEVANPEI-263-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170006480FD/output/BACKUP/Bamdir/170504_NS500823_HKG73BGX2_L001_V1D5HUMCNEVMMPEI-294-v14.sort.markdup.realign.recald.bam
	"""
	bams="""
	/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000558FD/output/BACKUP/Bamdir/170419_BS1001160019_CL100015943_L001_V1D5HUMRNETKOPEI-410-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000559FD/output/BACKUP/Bamdir/170419_BS1001160019_CL100015943_L001_V1D5HUMRNETKPPEI-411-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000568FD/output/BACKUP/Bamdir/170419_BS1001160019_CL100015943_L001_V1D5HUMRNETKXPEI-420-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000571FD/output/BACKUP/Bamdir/170419_BS1001160019_CL100015943_L001_V1D5HUMRNETKZPEI-422-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000601FD/output/BACKUP/Bamdir/170419_BS1001160019_CL100015943_L001_V1D5HUMRNETMCPEI-406-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000603FD/output/BACKUP/Bamdir/170421_BS1001160015_CL100015950_L001_V1D5HUMRNETMDPEI-402-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000605FD/output/BACKUP/Bamdir/170421_BS1001160015_CL100015950_L001_V1D5HUMRNETQCPEI-404-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000606FD/output/BACKUP/Bamdir/170421_BS1001160015_CL100015950_L001_V1D5HUMRNETQDPEI-405-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000608FD/output/BACKUP/Bamdir/170421_BS1001160015_CL100015950_L001_V1D5HUMRNETQFPEI-407-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000627FD/output/BACKUP/Bamdir/170421_BS1001160015_CL100015950_L001_V1D5HUMRNETQVPEI-424-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000648FD/output/BACKUP/Bamdir/170421_BS1001160015_CL100015950_L001_V1D5HUMRNETRJPEI-445-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000649FD/output/BACKUP/Bamdir/170421_BS1001160015_CL100015950_L001_V1D5HUMRNETRKPEI-446-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000668FD/output/BACKUP/Bamdir/170629_BS1001160241_CL100021692_L001_V1D5HUMRNFAUDPEI-413-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000671FD/output/BACKUP/Bamdir/170629_BS1001160241_CL100021378_L001_V1D5HUMRNFAUWPEI-434-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000674FD/output/BACKUP/Bamdir/170629_BS1001160241_CL100021378_L001_V1D5HUMRNFAVPPEI-405-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000676FD/output/BACKUP/Bamdir/170629_BS1001160241_CL100021378_L001_V1D5HUMRNFAVFPEI-443-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000677FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVKQPEI-419-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000680FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVKRPEI-420-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000691FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVKLPEI-414-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000695FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVKYPEI-427-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000702FD/output/BACKUP/Bamdir/170513_BS1001160015_CL100015948_L001_V1D5HUMRNEVKVPEI-424-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000704FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVLAPEI-429-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000748FD/output/BACKUP/Bamdir/170513_BS1001160015_CL100016039_L001_V1D5HUMRNEWDAPEI-426-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001023FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVKKPEI-413-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001042FD/output/BACKUP/Bamdir/170513_BS1001160015_CL100016039_L001_V1D5HUMRNEWBIPEI-430-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001058FD/output/BACKUP/Bamdir/170513_BS1001160015_CL100015948_L001_V1D5HUMRNEVZNPEI-431-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001060FD/output/BACKUP/Bamdir/170513_BS1001160015_CL100016039_L001_V1D5HUMRNEVKAPEI-403-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001063FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVIOPEI-413-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001064FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVJJPEI-434-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001067FD/output/BACKUP/Bamdir/170513_BS1001160015_CL100016039_L001_V1D5HUMRNEWCYPEI-424-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001070FD/output/BACKUP/Bamdir/170513_BS1001160015_CL100015948_L001_V1D5HUMRNEWBLPEI-433-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001074FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVIEPEI-403-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001089FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVIDPEI-402-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001105FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVJIPEI-433-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001153FD/output/BACKUP/Bamdir/170513_BS1001160015_CL100016039_L001_V1D5HUMRNEWCCPEI-402-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001169FD/output/BACKUP/Bamdir/170507_BS1001160019_CL100015960_L001_V1D5HUMRNEVJGPEI-431-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179001215FD/output/BACKUP/Bamdir/170629_BS1001160241_CL100021692_L001_V1D5HUMRNFATDPEI-435-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160009665FD/output/BACKUP/Bamdir/170221_NS500823_H53WJBGX2_L001_V1D5HUMCNEOMXPEI-211-v.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170005173FD/output/BACKUP/Bamdir/170422_TPNB500129_HKKKFBGX2_L001_HUMCNEULVPEI-263-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160007033FD/output/BACKUP/Bamdir/161104_TPNB500129_HV3T2BGXY_L001_V1D5HUMCAEFUYPEI-253-v.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160007559FD/output/BACKUP/Bamdir/161108_NS500823_H33L5BGX2_L001_V1D5HUMCAEGCAPEI-213-v.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160007942FD/output/BACKUP/Bamdir/161202_NS500823_H7NW7BGX2_L001_V1D5HUMCAEIXGPEI-252-v.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160007603FD/output/BACKUP/Bamdir/161212_NS500823_H7KNVBGX2_L001_V1D5HUMCAEJUUPEI-202-v13.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160009409FD/output/BACKUP/Bamdir/161220_TPNB500129_H7NHHBGX2_L001_V1D5HUMCAEKEEPEI-241-v106.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160004217FD/output/BACKUP/Bamdir/161224_TPNB500129_HF27LBGX2_L001_V1D5HUMCAEKNGPEI-239-v.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160007421FD/output/BACKUP/Bamdir/161231_TPNB500129_HF2HFBGX2_L001_V1D5HUMCAEKXNPEI-203-v106.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160009690FD/output/BACKUP/Bamdir/170103_TPNB500129_HCYW3BGX2_L001_V1D5HUMCAEKZXPEI-216-v106.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160007236FD/output/BACKUP/Bamdir/170109_TPNB500129_HCYNKBGX2_L001_V1D5HUMCAELNTPEI-259-v106.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160002650FD/output/BACKUP/Bamdir/170110_NS500823_HF2VWBGX2_L001_V1D5HUMCAELOLPEI-269-v.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160005892FD/output/BACKUP/Bamdir/170208_NS500823_H7YW3AFXX_L001_V1D5HUMCNENGQPEI-203-v106.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160009010FD/output/BACKUP/Bamdir/170211_NS500823_H5MLCBGX2_L001_V1D5HUMCNENPUPEI-219-v.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160002663FD/output/BACKUP/Bamdir/170218_NS500823_H53L2BGX2_L001_V1D5HUMCNEOJLPEI-224-v106.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170000362FD/output/BACKUP/Bamdir/170224_NS500823_H57CNBGX2_L001_V1D5HUMCNEOVWPEI-253-v.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160009403FD/output/BACKUP/Bamdir/170225_NS500823_H5MFYBGX2_L001_V1D5HUMCNEOYSPEI-295-v106.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160009149FD/output/BACKUP/Bamdir/170226_TPNB500129_H53F5BGX2_L001_V1D5HUMCNEOANPEI-228-v.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160005885FD/output/BACKUP/Bamdir/170312_TPNB500129_H5K5FBGX2_L001_HUMCNEQPPPEI-247-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/160002638FD/output/BACKUP/Bamdir/170321_TPNB500129_H5NHKBGX2_L001_HUMCNERLCPEI-227-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170002320FD/output/BACKUP/Bamdir/170325_TPNB500129_H5MTTBGX2_L001_V1D5HUMCNERWWPEI-261-v13.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170005652FD/output/BACKUP/Bamdir/170330_TPNB500129_H5MJJBGX2_L001_V1D5HUMCNESJDPEI-285-v13.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170005864FD/output/BACKUP/Bamdir/170403_TPNB500129_H5CKCBGX2_L001_HUMCNESUZPEI-221-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170005060FD/output/BACKUP/Bamdir/170408_TPNB500129_HKCHFBGX2_L001_HUMCNETBVPEI-287-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170006964FD/output/BACKUP/Bamdir/170410_TPNB500129_HKHCMBGX2_L001_HUMCNETEIPEI-202-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170005061FD/output/BACKUP/Bamdir/170411_TPNB500129_HKHJMBGX2_L001_V1D5HUMCNETGVPEI-220-v13.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170002604FD/output/BACKUP/Bamdir/170413_TPNB500129_HKHKJBGX2_L001_HUMCNETNEPEI-272-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170005063FD/output/BACKUP/Bamdir/170413_TPNB500129_HKHKJBGX2_L001_V1D5HUMCNETNDPEI-255-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170005400FD/output/BACKUP/Bamdir/170421_TPNB500129_HKHJKBGX2_L001_HUMCNEUIHPEI-293-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170006773FD/output/BACKUP/Bamdir/170421_TPNB500129_HKHJKBGX2_L001_HUMCNEUIGPEI-218-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170005156FD/output/BACKUP/Bamdir/170426_TPNB500129_HJV7HBGX2_L001_HUMCNEUVUPEI-295-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170006909FD/output/BACKUP/Bamdir/170428_TPNB500129_HVCNJBGX2_L001_HUMCNEVANPEI-263-v14.sort.markdup.realign.recald.bam
/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/170006480FD/output/BACKUP/Bamdir/170504_NS500823_HKG73BGX2_L001_V1D5HUMCNEVMMPEI-294-v14.sort.markdup.realign.recald.bam

	"""
	bams2="""
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/P-1_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/N3_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/L2-1_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/S2_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/P-5_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/P-7_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/P-10_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/L1-9_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/L2-6_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/L2-2_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/L2-9_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/N2_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/N4_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/N6_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/N5_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/L1-10_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/N_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/P-8_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/S3_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/P-6_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/P-3_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/S1_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/L1-1_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/S5_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/L1-7_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/L1-8_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/N8_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/P-2_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/P-9_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/N9_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/N1_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/L2-7_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/P-4_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/S4_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/N7_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/L1-6_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/L2-10_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/L1-4_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/L2-4_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/L2-5_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/L2-3_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/S6_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/L2-8_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/S7_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/L1-2_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/L1-3_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/L1-5_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/stand_sample/P_del_15bp.bam
	"""
	bam3="""
	/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000674FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000568FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001063FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000748FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000608FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160007603FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000558FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001023FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170005173FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160009149FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160002650FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170005156FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001153FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000691FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001074FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001215FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160005892FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001105FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170005652FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000676FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001089FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160009403FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001060FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170005864FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160009010FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170000362FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000648FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000704FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170005063FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000702FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160007559FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160007033FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000677FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000571FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001067FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170005061FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170006909FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160007236FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000559FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170005400FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170006964FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001042FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160007942FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000601FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000606FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160002638FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001064FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160007421FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001169FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170002604FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001070FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000695FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001058FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160004217FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000671FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170006773FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000668FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000627FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160009409FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160009690FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000649FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170006480FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170005060FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000680FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000605FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160009665FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160005885FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160002663FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170002320FD_del_15bp.abra.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000603FD_del_15bp.abra.bam
	"""
	bam4="""
	/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160009149FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170005864FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160009010FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170006909FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000627FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160005892FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000671FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160007603FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001063FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001215FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160005885FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001070FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000677FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000649FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160009690FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170005652FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000603FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001042FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160004217FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000691FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000674FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170005173FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000648FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001064FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000558FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170006480FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000676FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000608FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000695FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000559FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160007033FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160002663FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170005063FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001074FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170000362FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160007236FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000704FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160002650FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001169FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000748FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170005156FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001105FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000605FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000668FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170002320FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000571FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000568FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000601FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160002638FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170005060FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001067FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170006773FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160007942FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170002604FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160009403FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160009665FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001089FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001153FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170006964FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000702FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160009409FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160007559FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/160007421FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170005061FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001058FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/170005400FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000680FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179000606FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001023FD_del_15bp.bam
/mnt/NL200/zhouzhq/project/abra2_test/abra2/oncod_test/179001060FD_del_15bp.bam
	"""
	#/mnt/X500/farmers/chenzhx/AIO/analysis/2017/08/21/72samples/anadir/179000559FD/output/BACKUP/Bamdir/170419_BS1001160019_CL100015943_L001_V1D5HUMRNETKPPEI-411-v14.sort.markdup.realign.recald.bam
	#chr="chr7"          #1	26885310	26885311   7:55,242,454-55,242,489
	chr='7'          #7    55242464    55242479         7:55,242,454-55,242,489
	start=55242464
	end=55242479
	f_ref=open("deletion_ref_dug","w")
	f_mut=open("deletion_mut_dug","w")
	f_N_both=open("deletion_N_ref_and_mut","w")
	f_N_consistent=open("deletion_NoneConstent","w")
	#for bam in ["/mnt/NL200/wulj/work/CirculatingTumorDNA/abra2Test/output/170011503BD_170011503TD/cancer/4_realign_bam/170011503BD_170011503TD_cancer_sort_markdup_realignori.bam"]:
	for bam in bam0.split():
		result,mut_count,_=correct_under_del_rate(bam,chr,start,end,ref_genome)
		print "{}\t{}\t{}".format(bam.split("/")[12][:-1],result,mut_count)
		#print "{}\t{}\t{}".format(bam.split("/")[-1].split("_")[0][:-1],result,mut_count)
		#print "bam=",bam
		#print 'result=',result
		#print 'mut_count=',mut_count
	f_ref.close()
	f_mut.close()
	f_N_both.close()
	f_N_consistent.close()