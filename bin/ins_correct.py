#!/usr/bin/env python
# -*- coding:utf-8 -*-

from __future__ import division
import pysam,re,copy,sys,random,logging,multiprocessing
sys.dont_write_bytecode = True
debugging=False
f_ref,f_mut,f_N_both,f_N_consistent="","","",""
#HG="/mnt/NL200/prod/repos/BNC/BNC/program/NoahCare/db/alignment/tgp_phase2_flat/hs37d5.fa"

def mstch_count(a,b):
	logger=multiprocessing.get_logger()
	mismatch_count = 0
	a_len = len(a)
	b_len = len(b)
	if (a_len != b_len):
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

def handle_repeat_insert(chr,pos,insert_seq,HG):
	min_pos,max_pos=pos,pos
	i_len=len(insert_seq)
	l_start,l_end=pos,pos+i_len
	r_start,r_end=pos-i_len,pos
	left_loop_mark=True
	right_loop_mark=True
	while(left_loop_mark or right_loop_mark):
		if(left_loop_mark):
			del_seq_l=get_ref(chr,l_start-i_len,l_end-i_len,HG)
			for i in range(i_len-1,-1,-1):
				if(insert_seq[i]==del_seq_l[i]):
					min_pos-=1
				else:
					left_loop_mark=False
					break
				if(i==0):
					l_start-=i_len
					l_end-=i_len
		if(right_loop_mark):
			del_seq_r=get_ref(chr,r_start+i_len,r_end+i_len,HG)
			for i in range(i_len):
				if(insert_seq[i]==del_seq_r[i]):
					max_pos+=1
				else:
					right_loop_mark=False
					break
				if(i==(i_len-1)):
					r_start+=i_len
					r_end+=i_len
	return (min_pos,max_pos)

def modify_insert_type_read(cigar,read,del_parts,pos,blocks,del_blocks,chr,shift_start,shift_end,HG):
	logger=multiprocessing.get_logger()
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
		DM_len.append(int(item[:-1]))
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
		read_parts.append(read[cursor:sum[i]])
		cursor=sum[i]
	for i in range(len(D_M_array)):
		if("D" in D_M_array[i]):
			read_parts[i]=del_parts[0]
			if(len(del_parts[1:])>0):
				del_parts=del_parts[1:]
			else:
				break
	Insert_l_r_blocks=dict()
	M_match_blocks=dict()
	M_block_index=0
	for i in range(len(D_M_array)):
		if("M" in D_M_array[i]):
			M_match_blocks[i]=blocks[M_block_index]
			M_block_index+=1
	for i in range(len(D_M_array)):
		if("I" in D_M_array[i]):
			if(i>0 and i<(len(D_M_array)-1)):
				'''
				print 'D_M_array=',D_M_array
				print 'i=',i
				print 'len(M_match_blocks)=',len(M_match_blocks)
				print 'cigar=',cigar
				'''
				Insert_l_r_blocks[i]=[M_match_blocks[i-1],M_match_blocks[i+1]]
	Insert_delete_bool=dict()
	for i in Insert_l_r_blocks:
		pos_=Insert_l_r_blocks[i][0][1]
		if (pos_>=shift_start and pos_<=shift_end):
			if pos_!=pos:
				read_parts[i]=""
		else:
			read_parts[i]=""
	new_read_start,new_read_end=blocks[0][0],blocks[-1][1]
	if(M_index[0]>0):
		for i in range(M_index[0])[::-1]:
			if("S" in D_M_array[i]):
				if(len(read_parts[i])>10):
					S_len=DM_len[i]
					S_add_seq=get_ref(chr,new_read_start-S_len,new_read_start,HG)
					read_parts[i]=S_add_seq
				new_read_start-=len(read_parts[i])
	if(M_index[-1]!=len(D_M_array)-1):
		for i in range(M_index[-1]+1,len(D_M_array)):
			if("S" in D_M_array[i]):
				if(len(read_parts[i])>10):
					S_len=DM_len[i]
					S_add_seq=get_ref(chr,new_read_end,new_read_end+S_len,HG)
					read_parts[i]=S_add_seq
				new_read_end+=len(read_parts[i])
	new_read="".join(read_parts)
	return new_read,new_read_start,new_read_end


def get_modify_insert_D_parts(read,cigar,blocks,chr,pos,shift_start,shift_end,HG):
	is_mut=False
	pre_start,pre_end=blocks[0]
	del_parts=list()
	del_blocks=list()
	for i in range(1,len(blocks)):
		cur_start,cur_end=blocks[i]
		if(pre_end!=cur_start):
			del_blocks.append((pre_end,cur_start))
			del_seq=get_ref(chr,pre_end,cur_start,HG)
			del_parts.append(del_seq)
		pre_start,pre_end=cur_start,cur_end
	new_read,new_read_start,new_read_end=modify_insert_type_read(cigar,read,del_parts,pos,blocks,del_blocks,chr,shift_start,shift_end,HG)
	return new_read,new_read_start,new_read_end

def get_4_normal_insert_ref(chr,pos,shift_start,shift_end,insert_seq,HG):
	flk_len=400
	left_normal=get_ref(chr,shift_end-flk_len,shift_end,HG)
	right_normal=get_ref(chr,shift_start,shift_start+flk_len,HG)
	left_insert=get_ref(chr,pos-flk_len,pos,HG)+insert_seq+get_ref(chr,pos,shift_end,HG)
	right_insert=get_ref(chr,shift_start,pos,HG)+insert_seq+get_ref(chr,pos,pos+flk_len,HG)
	return left_normal,right_normal,left_insert,right_insert

def get_read_left_and_right(chr,read,shift_start,shift_end,R_start,R_end,insert_query_seq,HG):
	left_read,right_read="",""
	read_first,read_middle,read_three="","",""
	if(R_start<shift_start and R_end>shift_end):
		read_middle=read
		read_middle=read_middle[(shift_start-R_start):]
		read_middle=read_middle[:-(R_end-shift_end)]
		read_first=read[:(shift_start-R_start)]
		read_three=read[-(R_end-shift_end):]
		insert_len=len(insert_query_seq)
		left_len=shift_start-R_start
		right_len=R_end-shift_end
		left_read=read_first+read_middle
		right_read=read_middle+read_three
	if(R_end>=shift_start and R_end<=shift_end):
		left_remained=read[(shift_start-R_start):]
		left_read=read_first+left_remained+get_ref(chr,R_end+len(left_remained),shift_end,HG)
		right_read=left_remained
	if(R_start>=shift_start and R_start<=shift_end):
		right_remianed=read[:-(R_end-shift_end)]
		right_read=get_ref(chr,shift_start,R_end-len(right_remianed),HG)+right_remianed+read_three
		left_read=right_remianed
	return left_read,right_read

def is_insert(read_left,read_right,normal_left,normal_right,insert_normal_left,insert_normal_right,cigar):
	flag_other=0
	l_r_mismatch_count=None
	l_r_read_len=len(read_left)+len(read_right)
	l_r_mismatch_rate=None
	read_lest_len=len(read_left)
	read_right_len=len(read_right)
	same_len_normal_left=normal_left[-read_lest_len:]
	same_len_normal_right=normal_right[:read_right_len]
	same_len_insert_normal_left=insert_normal_left[-read_lest_len:]
	same_len_insert_normal_right=insert_normal_right[:read_right_len]
	left_normal_mismatch,left_insert_mismatch,right_normal_mismatch,right_insert_mismatch=0,0,0,0
	if(len(read_left)>0):
		left_normal_mismatch=mstch_count(read_left,same_len_normal_left)
		left_insert_mismatch=mstch_count(read_left,same_len_insert_normal_left)
	if(len(read_right)>0):
		right_normal_mismatch=mstch_count(read_right,same_len_normal_right)
		right_insert_mismatch=mstch_count(read_right,same_len_insert_normal_right)
	def print_varaible(f_handle):
		savedStdout = sys.stdout
		sys.stdout=f_handle
		print("cigar=",cigar)
		print("read_left=                  ",read_left)
		print("same_len_normal_left=       ",same_len_normal_left)
		print("same_len_insert_normal_left=",same_len_insert_normal_left)
		print("left_normal_mismatch=",left_normal_mismatch)
		print("left_insert_mismatch=",left_insert_mismatch)
		print("read_right=                  ",read_right)
		print("same_len_normal_right=       ",same_len_normal_right)
		print("same_len_insert_normal_right=",same_len_insert_normal_right)
		print("right_normal_mismatch=",right_normal_mismatch)
		print("right_insert_mismatch=",right_insert_mismatch)
		print(":"*60)
		sys.stdout = savedStdout
	print_varaible(sys.stdout) if debugging else False
	return_result=""
	f_out_handle=""
	if(left_normal_mismatch<left_insert_mismatch and right_normal_mismatch<right_insert_mismatch):
		l_r_mismatch_count=left_normal_mismatch+right_normal_mismatch
		return_result=0
		f_out_handle=f_ref
	elif(left_normal_mismatch>left_insert_mismatch and right_normal_mismatch>right_insert_mismatch):
		l_r_mismatch_count=left_insert_mismatch+right_insert_mismatch
		return_result=1
		f_out_handle=f_mut
	elif(left_normal_mismatch==left_insert_mismatch and right_normal_mismatch>right_insert_mismatch):
		l_r_mismatch_count=left_insert_mismatch+right_insert_mismatch
		return_result=1
		f_out_handle=f_mut
	elif(left_normal_mismatch==left_insert_mismatch and right_normal_mismatch<right_insert_mismatch):
		l_r_mismatch_count=left_normal_mismatch+right_normal_mismatch
		return_result=0
		f_out_handle=f_ref
	elif(left_normal_mismatch>left_insert_mismatch and right_normal_mismatch==right_insert_mismatch):
		l_r_mismatch_count=left_insert_mismatch+right_normal_mismatch
		return_result=1
		f_out_handle=f_mut
	elif(left_normal_mismatch<left_insert_mismatch and right_normal_mismatch==right_insert_mismatch):
		l_r_mismatch_count=left_normal_mismatch+right_normal_mismatch
		return_result=0
		f_out_handle=f_ref
	elif(left_normal_mismatch<left_insert_mismatch and right_normal_mismatch>right_insert_mismatch):
		if re.search(r"^\d+M$",cigar):
			l_r_mismatch_count=0
			return_result=0
			f_out_handle=f_ref
		else:
			if left_normal_mismatch+right_normal_mismatch < left_insert_mismatch+right_insert_mismatch:
				l_r_mismatch_count=left_normal_mismatch+right_normal_mismatch
				return_result=0
				f_out_handle=f_ref
			elif left_normal_mismatch+right_normal_mismatch > left_insert_mismatch+right_insert_mismatch:
				l_r_mismatch_count=left_insert_mismatch+right_insert_mismatch
				return_result=1
				f_out_handle=f_mut
			else:
				return_result=-2
				f_out_handle=f_N_consistent
				flag_other=1
		'''
		only_M_pat=re.compile(r'^[0-9]+M$')
		only_M_array=only_M_pat.findall(cigar)
		if(len(only_M_array)==1):
			return_result=0
			f_out_handle=f_ref
		return_result=-2
		f_out_handle=f_N_consistent
		'''
	elif(left_normal_mismatch>left_insert_mismatch and right_normal_mismatch<right_insert_mismatch):
		if re.search(r"^\d+M$",cigar):
			l_r_mismatch_count=0
			return_result=0
			f_out_handle=f_ref
		else:
			if left_normal_mismatch+right_normal_mismatch < left_insert_mismatch+right_insert_mismatch:
				l_r_mismatch_count=left_normal_mismatch+right_normal_mismatch
				return_result=0
				f_out_handle=f_ref
			elif left_normal_mismatch+right_normal_mismatch > left_insert_mismatch+right_insert_mismatch:
				l_r_mismatch_count=left_insert_mismatch+right_insert_mismatch
				return_result=1
				f_out_handle=f_mut
			else:
				return_result=-2
				f_out_handle=f_N_consistent
				flag_other=1
		'''
		only_M_pat=re.compile(r'^[0-9]+M$')
		only_M_array=only_M_pat.findall(cigar)
		if(len(only_M_array)==1):
			return_result=0
			f_out_handle=f_ref
		return_result=-2
		f_out_handle=f_N_consistent
		'''
	else:
		return_result=-1
		f_out_handle=f_N_both
		flag_other=1
	print_varaible(f_out_handle) if debugging else False
	if return_result==0 or return_result==1:
		l_r_mismatch_rate=(l_r_mismatch_count+0.0)/l_r_read_len
	return return_result,l_r_mismatch_rate,flag_other

def insert_should_consider_read(cigar,blocks,pos,shift_start,shift_end,R_start,R_end):
	if(R_end>=shift_start and R_end<=shift_end):
		return False
	if(R_start>=shift_start and R_start<=shift_end):
		return False
	M_start,M_end=blocks[0][0],blocks[-1][1]
	if(pos>=(M_start+1) and pos<=(M_end-1)):
		return True
	cigar_pat=re.compile(r'[0-9]+D|[0-9]+M|[0-9]+S|[0-9]+I')
	D_M_array=cigar_pat.findall(cigar)
	M_index=list()
	S_index=list()
	for i in range(len(D_M_array)):
		if("M" in D_M_array[i]):
			M_index.append(i)
		elif("S" in D_M_array[i]):
			S_index.append(i)
	DM_len=list()
	for item in D_M_array:
		DM_len.append(int(item[:-1]))
	if("S" in cigar):
		if(M_index[0]>0):
			for i in range(M_index[0])[::-1]:
				if("S" in D_M_array[i]):
					if(int(DM_len[i])<10):
						S_len=DM_len[i]
						M_start-=S_len
		if(M_index[-1]!=len(D_M_array)-1):
			for i in range(M_index[-1]+1,len(D_M_array)):
				if("S" in D_M_array[i]):
					if(int(DM_len[i])<10):
						S_len=DM_len[i]
						M_end+=S_len
	if(pos>=(M_start+1) and pos<=(M_end-1)):
		return True
	return False

def is_not_exact_I_len(cigar,insert_seq,shift_start,shift_end,blocks):
	cigar_pat=re.compile(r'[0-9]+D|[0-9]+M|[0-9]+S|[0-9]+I')
	D_M_array=cigar_pat.findall(cigar)
	M_index=list()
	I_index=list()
	for i in range(len(D_M_array)):
		if("M" in D_M_array[i]):
			M_index.append(i)
		elif("I" in D_M_array[i]):
			I_index.append(i)
	DM_len=list()
	for item in D_M_array:
		DM_len.append(int(item[:-1])) 
	pre_start,pre_end=blocks[0][0],blocks[0][1]
	I_count=0
	I_len=len(insert_seq)
	for i in range(1,len(blocks)):
		cur_start,cur_end=blocks[i]
		if(pre_end==cur_start):
			insert_len=DM_len[I_index[I_count]]
			if(insert_len!=I_len and cur_start>=shift_start and cur_start<=shift_end):
				return True
			I_count+=1
		pre_start,pre_end=cur_start,cur_end
	return False

def get_read_len(bam,chr,pos):
	samfile = pysam.AlignmentFile(bam, "rb")
	for read in samfile.fetch(chr,pos,pos+1):
		read_sequence=read.query_sequence
		read_len=len(read_sequence)
		samfile.close()
		return read_len

def insert_correct_main(chr,pos,insert_seq,bam,insert_ref_genome):
	ref_qname_mismatch_dict=dict()
	alt_qname_mismatch_dict=dict()
	ref_query_name_set=set()
	alt_query_name_set=set()
	HG=insert_ref_genome
	import pysam
	samfile = pysam.AlignmentFile(bam, "rb")
	loc="1	26885317	26885317"
	arr=loc.split("\t")
	read_len=get_read_len(bam,chr,pos)+4
	insert_count=0
	ref_count=0
	N_ref_and_mut=0
	N_nonconsistent=0
	result=""
	shift_start,shift_end=handle_repeat_insert(chr,pos,insert_seq,HG)
	left_normal,right_normal,left_insert_normal,right_insert_normal=get_4_normal_insert_ref(chr,pos,shift_start,shift_end,insert_seq,HG)
	other_read_count=0
	for read in samfile.fetch(chr,pos-read_len,pos+read_len):
		def print_varaible(f_handle):
			savedStdout = sys.stdout
			sys.stdout=f_handle
			print("blocks=",blocks)
			print("c_gar=",c_gar)
			print("shift_start=",shift_start)
			print("shift_end=",shift_end)
			print("pos=",pos)
			print("new_read=",new_read)
			print("R_start=",R_start)
			print("R_end=",R_end)
			print("left_normal=",left_normal)
			print("left_insert_normal=",left_insert_normal)
			print("right_normal=",right_normal)
			print("right_insert_normal=",right_insert_normal)
			print("left_read=",left_read)
			print("right_read=",right_read)
			print("insert_result_type=",insert_result_type)
			print("insert_count=",insert_count)
			print("ref_count=",ref_count)
			print("N_ref_and_mut=",N_ref_and_mut)
			print("N_nonconsistent=",N_nonconsistent)
			print("-"*60)
			sys.stdout = savedStdout
		if(read.is_unmapped or read.mate_is_unmapped):
			continue
		if(is_dup(read)):
			continue
		read_sequence=read.query_sequence
		r_query_name=read.query_name
		b_locks=read.get_blocks()
		c_gar=read.cigarstring
		blocks=read.get_blocks()
		if "M" not in c_gar:
			continue
		#if r_query_name=='NS500823:452:HTHKKBGX2:2:22209:6065:2414':
		#	print 'yes'
		#	print 'pos=',pos
		#	print 'c_gar=',c_gar
		#	print 'blocks=',blocks
		#	print 'read_sequence=',read_sequence
		m=re.search(r'^(\d+)M(\d+)I(\d+)M$',c_gar)
		if m:
			if blocks[0][1]==pos:
				if int(m.group(2))==len(insert_seq) and read_sequence[int(m.group(1)):int(m.group(1))+int(m.group(2))]:
					alt_query_name_set.add(r_query_name)
					alt_qname_mismatch_dict[r_query_name]=-1
				continue
			'''
			else:
				ref_query_name_set.add(r_query_name)
				ref_qname_mismatch_dict[r_query_name]=-1
			'''
		if re.search(r'(^\d+M\d+S$)|(^\d+S\d+M$)',c_gar) and (blocks[0][0]==pos or blocks[0][1]==pos ):
			m= re.search(r'^(\d+)S(\d+)M$',c_gar)
			if m:
				if read_sequence[:int(m.group(1))]==insert_seq[-int(m.group(1)):]:
					alt_query_name_set.add(r_query_name)
					alt_qname_mismatch_dict[r_query_name]=-1
					continue
			m= re.search(r'^(\d+)M(\d+)S$',c_gar)
			if m :
				if read_sequence[int(m.group(1)):]==insert_seq[:int(m.group(2))]:
					alt_query_name_set.add(r_query_name)
					alt_qname_mismatch_dict[r_query_name]=-1
					continue
		
		m= re.search(r'^(\d+)S\d+M\d+',c_gar)
		if m and blocks[0][0]==pos:
			if read_sequence[:int(m.group(1))]==insert_seq[-int(m.group(1)):]:
				alt_query_name_set.add(r_query_name)
				alt_qname_mismatch_dict[r_query_name]=-1
				continue
		m= re.search(r'\d+[A-Z]\d+M(\d+)S$',c_gar)
		if m and blocks[-1][1]==pos:
			if read_sequence[-int(m.group(1)):]==insert_seq[:int(m.group(1))]:
				alt_query_name_set.add(r_query_name)
				alt_qname_mismatch_dict[r_query_name]=-1
				continue
		
		m=re.search(r'^\d+M$',c_gar)
		if m and blocks[0][0] <pos and blocks[0][1]>pos:
			ref_query_name_set.add(r_query_name)
			ref_qname_mismatch_dict[r_query_name]=-1
			continue
		b_start=blocks[0][0]
		b_end=blocks[-1][1]
		if not (b_start<pos and pos<b_end):
			continue
		cigar_pat=re.compile(r'[0-9]+M')
		D_array=cigar_pat.findall(c_gar)
		if "I" in c_gar and len(D_array)<2:
			other_read_count+=1
			continue
		D_len=int(D_array[0][:-1])
		
		if("I" in c_gar and len(blocks)>=2):
			if(is_not_exact_I_len(c_gar,insert_seq,shift_start,shift_end,blocks)):
				other_read_count+=1
				continue
		block_left_start,block_right_end=blocks[0][0],blocks[-1][1]
		other_insert_pos=block_left_start+D_len
		min_pos=min(pos,other_insert_pos)
		max_pos=max(pos,other_insert_pos)
		min_max_10bp=get_ref(chr,min_pos-10,max_pos+10,HG)
		new_read,R_start,R_end=get_modify_insert_D_parts(read_sequence,c_gar,blocks,chr,pos,shift_start,shift_end,HG)
		if(R_start>=shift_start and R_end<=shift_end):
			continue
		should_conside=insert_should_consider_read(c_gar,blocks,pos,shift_start,shift_end,R_start,R_end)
		if(not should_conside):
			other_read_count+=1
			continue
		left_read,right_read=get_read_left_and_right(chr,new_read,shift_start,shift_end,R_start,R_end,insert_seq,HG)
		insert_result_type,mismatch_rate,flag_other=is_insert(left_read,right_read,left_normal,right_normal,left_insert_normal,right_insert_normal,c_gar)
		other_read_count+=flag_other
		f_out_handle=""
		if(insert_result_type==1): #insert_result_type
			alt_qname_mismatch_dict[r_query_name]=mismatch_rate
			insert_count+=1
			alt_query_name_set.add(r_query_name)
			f_out_handle=f_mut
		elif(insert_result_type==0):
			ref_qname_mismatch_dict[r_query_name]=mismatch_rate
			ref_count+=1
			ref_query_name_set.add(r_query_name)
			f_out_handle=f_ref
		elif(insert_result_type==-1):
			N_ref_and_mut+=1
			f_out_handle=f_N_both
		elif(insert_result_type==-2):
			N_nonconsistent+=1
			f_out_handle=f_N_consistent
		print_varaible(f_out_handle) if debugging else False
		print_varaible(sys.stdout) if debugging else False
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
	if(len(alt_query_name_set)+len(ref_query_name_set)+other_read_count!=0):
		#result=insert_count/(insert_count+ref_count)
		result=len(alt_query_name_set)/(len(alt_query_name_set)+len(ref_query_name_set)+other_read_count)
	else:
		result=0
	if(debugging):
		print("insert_count/(insert_count+ref_count)=",result)
	samfile.close()
	#for i in alt_query_name_set:
	#	print 	i 
	return result,len(alt_query_name_set),len(alt_query_name_set)+len(ref_query_name_set)+other_read_count

if __name__=="__main__":
	debugging=False
	f_ref=open("insertion_ref_dug","w")
	f_mut=open("insertion_mut_dug","w")
	f_N_both=open("insertion_N_ref_and_mut","w")
	f_N_consistent=open("insertion_NoneConstent","w")
	'''
	bam="/mnt/NL200/prod/workspace/IFA20161021001/OncoD/output/160005583BCD_160005583BPD/cancer/5_recal_bam/160005583BCD_160005583BPD_cancer_sort_markdup_realign_recal.bam"
	insert_ref_genome="/mnt/NL200/prod/repos/BNC/BNC/program/NoahCare/db/alignment/tgp_phase2_flat/hs37d5.fa"
	loc="10	89726688	89726688"
	insert_seq="T"
	'''
	bam='/mnt/NL200_1/prod/workspace/IFA20170529004/OncoD/output/170007656B2CD_170007656TD/cancer/5_recal_bam/170007656B2CD_170007656TD_cancer_sort_markdup_realign_recal_ori.bam'
	insert_ref_genome="/mnt/NL200/prod/repos/BNC/BNC/program/NoahCare/db/alignment/tgp_phase2_flat/hs37d5.fa"
	loc="1	204521408	204521408"
	insert_seq="A"
	arr=loc.split("\t")
	chr=str(arr[0])
	pos=int(arr[1])
	insert_correct_main(chr,pos,insert_seq,bam,insert_ref_genome)
	f_ref.close()
	f_mut.close()
	f_N_both.close()
	f_N_consistent.close()