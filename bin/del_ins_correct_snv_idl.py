#!/usr/bin/env python
# -*- coding:utf-8 -*-

import csv,os,argparse,string,sys,multiprocessing,logging,functools,time
sys.dont_write_bytecode = True
csv.register_dialect("line_terminator",lineterminator="\n")
VERSION='1.0.0'
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
	line=list(line)
	#if(reader.line_num==1):
	#	writer.writerow(line)
	#	continue
	Chr_index,Start_index,End_index,Ref_index,Alt_index,MutType_index,Case_var_freq_index=index_list
	chr=line[Chr_index]
	start=int(line[Start_index])
	end=int(line[End_index])
	mode=line[MutType_index]
	del_seq=line[Ref_index]
	insert_seq=line[Alt_index]
	result=None
	o_result=line[Case_var_freq_index]
	n_result=None
	count_line=",".join(map(str,[chr,start,end,mode,del_seq,insert_seq]))
	mut_count=0
	all_count=None
	if(mode=='del'):
		from del_correct import correct_under_del_rate as correct_under_del_rate
		n_result,mut_count,all_count=correct_under_del_rate(bam,chr,start,end,ref_genome)
	elif(mode=='ins'):
		from ins_correct import insert_correct_main as insert_correct
		n_result,mut_count,all_count=insert_correct(chr,start,insert_seq,bam,ref_genome)
	if mode=='del' or mode=='ins':
		if n_result > 0:
			if float(n_result)*100>float(o_result):
				result="00"+str(n_result*100)
			elif float(n_result)*100<float(o_result):
				result="0"+str(n_result*100)
			else:
				result=str(n_result*100)
			line[Case_var_freq_index]=result
		else:
			line[Case_var_freq_index]=o_result
	if((mode=='del' or mode=='ins') and result is not None):
		count_line=count_line+","+str(mut_count)+","+str(result)
		return (line,count_line)
	else:
		return (line,None)

def del_ins_correct_main(csv_file,bam,ref_genome,f_out_name,out_count_file,pro_num):
	col_names=['Chr','Start','End','Ref','Alt','MutType','Case_var_freq']
	f_count_h=open(out_count_file,"w")
	f_count_header=",".join(["chr","start","end","model","ref","alt","Crect_InDel_reds_count","Crect_InDel_reds_percent"])
	print>>f_count_h,f_count_header
	csvfile = file(csv_file, 'rb') 
	out_csvfile=file(f_out_name,'wb')
	reader = csv.reader(csvfile)
	writer = csv.writer(out_csvfile,dialect='line_terminator')
	header_line=next(reader)
	#Chr_index,Start_index,End_index,Ref_index,Alt_index,MutType_index,Case_var_freq_index=find_col_index(header_line,col_names)
	index_list=find_col_index(header_line,col_names)
	writer.writerow(list(header_line))
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