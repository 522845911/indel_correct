#!/usr/bin/env python
# -*- coding:utf-8 -*-

import csv,os,argparse,string,BedAnnoVar,sys,multiprocessing,logging,functools,time,signal, re
sys.dont_write_bytecode = True
csv.register_dialect("line_terminator",lineterminator="\n")
VERSION='1.0.1'
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

def find_col_index(header,col_names_dict):
	arr=list(header)
	index_list=list()
	for i in col_names_dict: 
		for j in arr:
			if i==j.strip("#"):
				index_list.append(arr.index(j))
	return index_list
@LogExc
def process_func(line,index_list,bam,ref_genome):
	chr_list=map(str,list(range(1,23)))+['X','Y'] 
	line=list(line)
	chr_index,pos_index,ref_index,alt_index,AF_index,altdep_index,totaldep_index,caseMuTect_index=index_list
	#if(reader.line_num==1):
	#	writer.writerow(line)
	#	continue
	chr=line[chr_index]
	pos=line[pos_index]
	ref=line[ref_index]
	alt=line[alt_index]
	caseMuTect=line[caseMuTect_index]
	indel_len = abs(len(ref)-len(alt))
	std_var=None
	try:
		std_var=BedAnnoVar.BedAnnoVar(chr,pos,ref,alt)
	except Exception,e:
		try:
			multiprocessing.get_logger().exception("<<<<<the *hot_spots.tsv file content has some probles!,pls check it!>>>>>")
			os.kill(os.getppid(), signal.SIGKILL)
		except OSError, e:
			pass
		sys.exit(1)
	
	start=int(std_var.pos)
	end=int(std_var.end)
	mode=std_var.imp
	if re.search(r'\|CX$',caseMuTect):
		mode="MuTect"
	if indel_len > 35:
		mode = "discard"
	del_seq="."
	insert_seq="."
	std_ref=std_var.ref
	if mode=="ins":
		insert_seq=std_var.alt
	elif mode=="del":
		del_seq=std_var.ref
	o_result=line[AF_index]
	count_line="\t".join(map(str,[chr,start,end,mode,del_seq,insert_seq]))
	mut_count="."
	all_count="."
	n_result=None
	result="." 
	if(chr in chr_list):
		if(mode=='del'):
			from del_correct import correct_under_del_rate as correct_under_del_rate
			n_result,mut_count,all_count=correct_under_del_rate(bam,chr,start,end,ref_genome)
		elif(mode=='ins'):
			from ins_correct import insert_correct_main as insert_correct
			n_result,mut_count,all_count=insert_correct(chr,start,insert_seq,bam,ref_genome)
		if mode=='del' or mode=='ins':
			if n_result > 0:
				if float(n_result)>float(o_result):
					result="00"+str(n_result)
				elif float(n_result)<float(o_result):
					result="0"+str(n_result)
				else:
					result=str(n_result)
	line.append(result)
	line.append(mut_count)
	line.append(all_count)
	if((mode=='del' or mode=='ins') and result is not None):
		count_line=count_line+"\t"+str(mut_count)+"\t"+str(result)
		return (line,count_line)
	else:
		return (line,None)

def del_ins_correct_main(csv_file,bam,ref_genome,f_out_name,out_count_file,pro_num):
	#col_names=['CHR','POS','REF','ALT','caseAltAF','caseAltDep','caseDep', 'caseMuTect']
	col_names=['CHR','POS','REF','ALT','caseAltAF','caseAltDep','caseRawDep', 'caseMuTect']
	f_count_h=open(out_count_file,"w")
	f_count_header="\t".join(["chr","start","end","model","ref","alt","Crect_InDel_reds_count","Crect_InDel_reds_percent"])
	print>>f_count_h,f_count_header
	csvfile = file(csv_file, 'rb')
	out_csvfile=file(f_out_name,'wb')
	reader = csv.reader(csvfile,delimiter="\t")
	writer = csv.writer(out_csvfile,delimiter="\t",dialect='line_terminator')
	header_line=next(reader)
	writer.writerow(list(header_line)+["caseCorrAF", "caseCorrAltDep", "caseCorrDep"])
	#chr_index,pos_index,ref_index,alt_index,AF_index,altdep_index,totaldep_index=find_col_index(header_line,col_names)
	index_list=find_col_index(header_line,col_names)
	pool = multiprocessing.Pool(processes=pro_num)
	result = list()
	for line in reader:
		result.append(pool.apply_async(process_func, (line, index_list,bam,ref_genome )))
	pool.close()
	pool.join()
	for res in result:
		(csv_line,csv_count_line)=res.get()
		writer.writerow(csv_line)
		if csv_count_line is not None:
			print>>f_count_h,csv_count_line
	csvfile.close()
	out_csvfile.close()
	f_count_h.close()

if __name__=="__main__":
	parser = argparse.ArgumentParser()  
	parser.add_argument('-c','--csv', help='the input csv file',required=True)
	parser.add_argument('-b', '--bam', help='the input bam file',required=True)
	parser.add_argument('-r', '--ref', help='the input reference genome',required=True)
	parser.add_argument('-o', '--out', help='the output csv file',required=True)
	parser.add_argument('-i', '--out_inDel_reads_count', help='the output in_del support reads count file',required=True)
	parser.add_argument('-p', '--pro_num', help='the process num',required=True,type=int)
	parser.add_argument('-l', '--logging_file', help='the logging file',required=True)
	args = parser.parse_args()
	csv_file=os.path.abspath(args.csv) 
	bam_file=os.path.abspath(args.bam)
	ref_file=os.path.abspath(args.ref)
	out_file=os.path.abspath(args.out)
	out_count_file=os.path.abspath(args.out_inDel_reads_count)
	pro_num=args.pro_num
	log_file=os.path.abspath(args.logging_file)
	f=logging.Formatter('[%(levelname)s %(processName)s %(asctime)s %(funcName)s] %(message)s')
	h=logging.FileHandler(log_file, 'w')
	h.setFormatter(f)
	logger=multiprocessing.get_logger()
	logger.addHandler(h)
	logger.setLevel(logging.INFO)
	logger.info("The program Version is %s",VERSION)
	logger.info('\nprogram del_ins_correct_nc_hot started in %s with command: %s', os.getcwd(), ' '.join(sys.argv))
	logger.info('\n[ the number of process is  %s]',pro_num)
	logger.info('\n[ start time: %s]',time.asctime( time.localtime(time.time()) ))
	del_ins_correct_main(csv_file,bam_file,ref_file,out_file,out_count_file,pro_num)
	logger.info('\n[ end time: %s]',time.asctime( time.localtime(time.time()) ))
